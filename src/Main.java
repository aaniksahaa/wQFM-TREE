package src;

import java.io.FileWriter;
import java.io.IOException;

import src.InitialPartition.ConsensusTreePartition;
import src.InitialPartition.IMakePartition;
import src.PreProcessing.GeneTrees;

import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;

/**
 * Main entry point for wQFM-TREE algorithm implementation.
 * 
 * wQFM-TREE is a quartet-based species tree inference method that applies the weighted
 * Quartet Fiduccia-Mattheyses (wQFM) algorithm directly to gene trees without requiring
 * decomposition into induced quartets. This implementation follows the methodology
 * described in the research paper "wQFM-TREE: highly accurate and scalable quartet-based
 * species tree inference from gene trees".
 * 
 * The algorithm uses a divide-and-conquer approach with two key innovations:
 * 1. Gene tree consensus-based initial bipartition computation (Algorithm 1 in paper)
 * 2. Direct scoring of candidate bipartitions from gene trees (Algorithm 2 in paper)
 * 
 * Time complexity: O(n³k log n) under balanced partitioning assumptions
 * where n = number of taxa, k = number of gene trees
 */
// 5:18
public class Main {

    /**
     * Main method that orchestrates the complete wQFM-TREE algorithm execution.
     * 
     * Command line arguments:
     * args[0] - Input file path containing gene trees in Newick format
     * args[1] - Consensus tree file path for initial bipartition construction
     * args[2] - Output file path for the resulting species tree
     * args[3] - Non-quartet type ("A" or "B") determining unresolved quartet handling
     * 
     * The non-quartet type corresponds to different formulations for computing
     * w(U^[g]) (sum of weights of unresolved quartets) as described in Section 2.5.3
     * of the paper.
     */
    public static void main(String[] args) throws IOException {

        // Validate command line arguments
        if(args.length < 4){
            System.out.println("Specify all file paths and non quartet type");
            System.exit(-1);
        }
        String inputFilePath = args[0];
        String consensusFilePath = args[1];
        String outputFilePath = args[2];

        String nonQuartetType = args[3];
        
        // Configure non-quartet handling strategy
        // This affects how unresolved quartets in polytomy nodes are processed
        // as described in Section 2.5.3 of the paper
        if(nonQuartetType.equals("A")){
            Config.NON_QUARTET_TYPE = Config.NonQuartetType.A;
        }else if(nonQuartetType.equals("B")){
            Config.NON_QUARTET_TYPE = Config.NonQuartetType.B;
        }else{
            System.out.println("Specify non quartet type as A or B");
            System.exit(-1);
        }

        // Initialize timing measurements for performance analysis
        ThreadMXBean threadMXBean = ManagementFactory.getThreadMXBean();
        long time_1 = System.currentTimeMillis(); //calculate starting time
        long cpuTimeBefore = threadMXBean.getCurrentThreadCpuTime();

        /**
         * STEP 1: Gene Tree Preprocessing
         * Read and parse input gene trees, building internal data structures.
         * This corresponds to the input processing phase mentioned in Section 2.
         */
        GeneTrees trees = new GeneTrees(inputFilePath);
        var taxaMap = trees.readTaxaNames();

        /**
         * STEP 2: Initialize Consensus Tree Partition Maker
         * Creates the initial bipartition constructor using consensus tree approach.
         * This implements Algorithm 1 from the paper: "Creation of initial bipartition
         * for a taxa set which will be further refined by the FM algorithm"
         * 
         * The consensus tree is used to generate candidate initial bipartitions
         * for each divide step, as described in Section 2.4.
         */
        ConsensusTreePartition consensusTreePartition = new ConsensusTreePartition(consensusFilePath, taxaMap, trees);
        
        // Complete gene tree reading with distance matrix for consensus-based partitioning
        trees.readGeneTrees(consensusTreePartition.dist);

        // Set the partition maker interface for the main algorithm
        IMakePartition  partitionMaker = consensusTreePartition;

        /**
         * STEP 3: Execute wQFM-TREE Algorithm
         * Initialize and run the main wQFM-TREE algorithm (QFM2 class).
         * This implements the complete divide-and-conquer framework described
         * in the paper, combining:
         * - Initial bipartition generation using consensus trees
         * - Iterative bipartition refinement using FM algorithm
         * - Direct scoring from gene trees without quartet enumeration
         */
        var qfm = new QFM2(trees, trees.taxa, partitionMaker);
        
        // Execute the main algorithm - returns the inferred species tree
        var spTree = qfm.runWQFM();

        /**
         * STEP 4: Output Results
         * Write the resulting species tree in Newick format to the specified output file.
         */
        FileWriter writer = new FileWriter(outputFilePath);
        writer.write(spTree.getNewickFormat());
        writer.close();
        
        // Performance reporting
        long cpuTimeAfter = threadMXBean.getCurrentThreadCpuTime();

        long time_del = System.currentTimeMillis() - time_1;
        long minutes = (time_del / 1000) / 60;
        long seconds = (time_del / 1000) % 60;
        System.out.format("\nElapsed Time taken = %d ms ==> %d minutes and %d seconds.\n", time_del, minutes, seconds);

        long cpuTimeUsed = cpuTimeAfter - cpuTimeBefore;

        seconds = cpuTimeUsed / 1_000_000_000;
        minutes = seconds / 60;
        seconds = seconds % 60;

        System.out.println("CPU time used: " + minutes + " minutes, " + seconds + " seconds");
    }
}
