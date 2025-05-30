package src;

import java.util.ArrayList;

import src.DSPerLevel.BookKeepingPerLevel;
import src.DSPerLevel.BookKeepingPerLevelv2;
import src.DSPerLevel.TaxaPerLevelWithPartition;
import src.InitialPartition.IMakePartition;
import src.PreProcessing.GeneTrees;
import src.Taxon.DummyTaxon;
import src.Taxon.RealTaxon;
import src.Tree.Tree;
import src.Tree.TreeNode;

/**
 * QFM2: Core implementation of the wQFM-TREE algorithm.
 * 
 * This class implements the weighted Quartet Fiduccia-Mattheyses algorithm applied
 * directly to gene trees, as described in the research paper. The algorithm follows
 * a divide-and-conquer approach where each divide step creates a bipartition of the
 * taxa set, recursively solves subproblems, and combines the results.
 * 
 * Key innovations implemented here:
 * 1. Direct application to gene trees without quartet enumeration
 * 2. Consensus tree-based initial bipartition generation (Algorithm 1)
 * 3. Direct bipartition scoring from gene trees (Algorithm 2)
 * 4. FM algorithm for iterative bipartition refinement
 * 
 * The algorithm maintains dummy taxa throughout the divide-and-conquer process
 * to enable proper combination of subproblem solutions, as described in Section 2.3
 * of the paper regarding "Tree structure of a dummy taxon and weighting scheme".
 */
public class QFM2 {
    
    // Core algorithm components
    public RealTaxon[] realTaxa;           // Original taxa from input gene trees
    public IMakePartition initPartition;   // Initial bipartition generator (Algorithm 1)
    public GeneTrees geneTrees;            // Input gene trees
    private int level;                     // Current recursion depth in divide-and-conquer

    // Numerical precision threshold for gain comparisons in FM algorithm
    static double EPS = 1e-5;

    /**
     * Constructor for wQFM-TREE algorithm.
     * 
     * @param trees Input gene trees in processed format
     * @param realTaxa Array of real taxa (leaves in gene trees)
     * @param initPartition Strategy for generating initial bipartitions (Algorithm 1)
     */
    public QFM2(GeneTrees trees, RealTaxon[] realTaxa, IMakePartition initPartition){
        this.realTaxa = realTaxa;
        this.initPartition = initPartition;
        this.geneTrees = trees;
    }

    /**
     * Main entry point for wQFM-TREE algorithm execution.
     * 
     * Initiates the divide-and-conquer process starting from the complete taxa set.
     * This corresponds to the top-level call in the recursive framework described
     * in the paper.
     * 
     * @return Complete species tree inferred from input gene trees
     */
    public Tree runWQFM(){
        this.level = 0;

        // Create initial bipartition for the complete taxa set
        // No dummy taxa exist at the root level
        var y = initPartition.makePartition(realTaxa, new DummyTaxon[0], true);
        var x = new TaxaPerLevelWithPartition(realTaxa, new DummyTaxon[0], y.realTaxonPartition, y.dummyTaxonPartition, realTaxa.length);
        
        // Initialize bookkeeping structure for scoring calculations
        // This handles the direct scoring from gene trees (Algorithm 2)
        BookKeepingPerLevelv2 initialBook = new BookKeepingPerLevelv2(geneTrees);
        initialBook.resetBookKeeping(x);

        return recurse(initialBook);
    }

    /**
     * Recursive divide-and-conquer method implementing the core wQFM-TREE algorithm.
     * 
     * Each call to this method handles one "divide step" as described in the paper:
     * 1. Generate initial bipartition using consensus tree (Algorithm 1)
     * 2. Refine bipartition using iterative FM algorithm
     * 3. Create two subproblems with dummy taxa
     * 4. Recursively solve subproblems
     * 5. Combine results by merging trees at dummy taxon positions
     * 
     * This implements the divide-and-conquer framework described in Section 2.5
     * and follows the tree structure methodology from Section 2.3.
     * 
     * @param book BookKeeping structure containing current taxa set and partition state
     * @return Rooted tree for the current subproblem
     */
    private Tree recurse(
        BookKeepingPerLevelv2 book
    ){

        this.level++;
        int itrCount = 0;

        // Singleton partition control based on recursion depth
        // This prevents excessive fragmentation in deep recursion levels
        boolean allowSingleton = Config.ALLOW_SINGLETON;

        if(allowSingleton){
            if(Config.USE_LEVEL_BASED_SINGLETON_THRESHOLD){
                // Level-based threshold: prevent singletons after certain recursion depth
                if(this.level > Config.MAX_LEVEL_MULTIPLIER * this.realTaxa.length){
                    System.out.println("Made false");
                    System.out.println("rts : " + book.taxas.realTaxonCount + " dts : " + book.taxas.dummyTaxonCount);
                    System.out.println("part[0]: " + book.taxas.getTaxonCountInPartition(0) + " part[1]: " + book.taxas.getTaxonCountInPartition(1));
                    allowSingleton = false;
                }
            }
            else{
                // Dummy taxa nesting depth threshold
                // Prevents excessive nesting of dummy taxa structures
                int maxDepth = 0;
                for(var x : book.taxas.dummyTaxa){
                    maxDepth = Math.max(maxDepth, x.nestedLevel);
                }
                if(maxDepth > Config.SINGLETON_THRESHOLD * this.realTaxa.length){
                    allowSingleton = false;
                }
            }
        }
        
        /**
         * ITERATIVE BIPARTITION REFINEMENT USING FM ALGORITHM
         * 
         * This loop implements the iterative improvement strategy described in the paper.
         * Each iteration attempts to improve the current bipartition by transferring
         * taxa between partitions. The loop continues until no improving moves are found.
         * 
         * The FM algorithm is a key component that distinguishes wQFM from other methods
         * like QMC, requiring specialized adaptation for direct application to gene trees.
         */
        while(oneInteration(book) ){
            itrCount++;
            if(itrCount > Config.MAX_ITERATION){
                System.out.println("Max iteration reached");
                break;
            }
        }
        System.out.println("level : " + level);
        System.out.println( "#iterations: " + itrCount);

        // Prepare for divide step: create two subproblems
        Tree[] trees = new Tree[2];

        /**
         * DIVIDE STEP: Create two subproblems with dummy taxa
         * 
         * This corresponds to the divide step described in Section 2 of the paper.
         * Each subproblem gets a unique dummy taxon representing the taxa in the
         * opposite partition. This enables proper combination in the conquer step.
         */
        var x = book.divide(initPartition, allowSingleton);
        int i = 0;
        int[] dummyIds = new int[2];
        
        // Recursively solve each subproblem
        for(var taxaWPart : x){
            book.resetBookKeeping(taxaWPart);
            if(book.taxas.smallestUnit){
                // Base case: create star tree for small subproblems
                trees[i] = book.taxas.createStar();
            }
            else{
                // Recursive case: further divide the subproblem
                trees[i] = recurse(book);
            }
            // Track dummy taxon IDs for tree combination
            dummyIds[i++] = taxaWPart.dummyTaxa[taxaWPart.dummyTaxonCount - 1].id;
        }

        this.level--;

        /**
         * CONQUER STEP: Combine subproblem solutions
         * 
         * This implements the tree combination strategy described in the paper.
         * The two subtrees are merged by connecting them at the positions of
         * their respective dummy taxa, effectively replacing dummy taxa with
         * the actual subtrees they represent.
         */
        TreeNode[] dtNodes = new TreeNode[2];

        // Locate dummy taxon nodes in both subtrees
        for(i = 0; i < 2; ++i){
            for(var node : trees[i].nodes){
                if(node.info.dummyTaxonId == dummyIds[i]){
                    dtNodes[i] = node;
                }
            }
            if(dtNodes[i] == null){
                System.out.println("Error: Dummy node not found");
                System.exit(-1);
            }
            if(dtNodes[i].childs != null){
                System.out.println("Error: Dummy node should be leaf");
                System.exit(-1);
            }
        }
        
        // Perform tree combination by re-rooting and grafting
        trees[0].reRootTree(dtNodes[0]);
        if(dtNodes[0].childs.size() > 1){
            System.out.println("Error: Dummy node should have only one child after reroot");
            System.exit(-1);
        }
        
        // Graft the first subtree onto the second at the dummy taxon position
        dtNodes[0].childs.get(0).setParent(dtNodes[1].parent);
        dtNodes[1].parent.childs.remove(dtNodes[1]);
        dtNodes[1].parent.childs.add(dtNodes[0].childs.get(0));
        trees[1].nodes.addAll(trees[0].nodes);

        return trees[1];
    }

    /**
     * Helper class for tracking taxon swaps during FM algorithm iterations.
     * 
     * Each swap represents the transfer of a taxon from one partition to another,
     * along with the associated gain in bipartition score.
     */
    static class Swap{
        public int index;        // Index of the taxon being swapped
        public boolean isDummy;  // Whether the taxon is a dummy taxon
        public double gain;      // Score improvement from this swap

        public Swap(int i, boolean id, double g){
            this.index = i;
            this.isDummy = id;
            this.gain = g;
        }
    }
    
    /**
     * Selects and executes the best taxon swap for FM algorithm iteration.
     * 
     * This method implements the core of the FM algorithm by:
     * 1. Evaluating all possible taxon transfers
     * 2. Selecting the transfer with maximum gain
     * 3. Executing the transfer and locking the moved taxon
     * 
     * The gain calculation uses the direct scoring method from Algorithm 2,
     * avoiding the need to enumerate quartets explicitly.
     * 
     * @param book Current bookkeeping state
     * @param rtGains Gains for real taxon transfers [taxon][partition]
     * @param dtGains Gains for dummy taxon transfers [taxon]
     * @param rtLocked Lock status for real taxa (prevents re-swapping)
     * @param dtLocked Lock status for dummy taxa
     * @return Swap object describing the executed transfer, or null if no valid swap
     */
    public static Swap swapMax(BookKeepingPerLevelv2 book, double[][] rtGains, double[] dtGains, boolean[] rtLocked, boolean[] dtLocked){

        int maxGainIndex = -1;
        double maxGain = 0;

        // Evaluate real taxon transfers
        for(int i = 0; i < book.taxas.realTaxonCount; ++i){
            if(rtLocked[i]) continue;
            int partition = book.taxas.inWhichPartitionRealTaxonByIndex(i);
            
            // Ensure partition size constraints are maintained
            if(book.taxas.getTaxonCountInPartition(partition) > 2 || (book.allowSingleton && book.taxas.getTaxonCountInPartition(partition) > 1) ){
                if(maxGainIndex == -1){
                    maxGain = rtGains[i][partition];
                    maxGainIndex = i;
                }
                else if(maxGain < rtGains[i][partition]){
                    maxGain = rtGains[i][partition];
                    maxGainIndex = i;
                }
            }
        }

        boolean dummyChosen = false;

        // Evaluate dummy taxon transfers
        for(int i = 0; i < book.taxas.dummyTaxonCount; ++i){
            if(dtLocked[i]) continue;
            int partition = book.taxas.inWhichPartitionDummyTaxonByIndex(i);
            
            // Ensure partition size constraints are maintained
            if(book.taxas.getTaxonCountInPartition(partition) > 2 || (book.allowSingleton && book.taxas.getTaxonCountInPartition(partition) > 1)){
                if(maxGainIndex == -1){
                    maxGain = dtGains[i];
                    maxGainIndex = i;
                    dummyChosen = true;
                }
                else if(maxGain < dtGains[i]){
                    maxGain = dtGains[i];
                    maxGainIndex = i;
                    dummyChosen = true;
                }
            }
        }

        // Return null if no valid swap found
        if(maxGainIndex == -1) return null;

        // Execute the selected swap and lock the taxon
        book.swapTaxon(maxGainIndex, dummyChosen);
        if(dummyChosen){
            dtLocked[maxGainIndex] = true;
        }
        else{
            rtLocked[maxGainIndex] = true;
        }

        return new Swap(maxGainIndex, dummyChosen, maxGain);
    }

    /**
     * Performs one iteration of the Fiduccia-Mattheyses (FM) algorithm for bipartition refinement.
     * 
     * This method implements the iterative improvement strategy described in Section 2 of the paper.
     * The FM algorithm is a key component that distinguishes wQFM from other quartet-based methods
     * like QMC. Each iteration:
     * 
     * 1. Evaluates all possible taxon transfers between partitions
     * 2. Selects and executes the transfer with maximum gain
     * 3. Locks the transferred taxon to prevent cycling
     * 4. Repeats until no more beneficial transfers exist
     * 5. Determines the optimal stopping point based on cumulative gain
     * 
     * The gain calculations are performed using the direct scoring method from Algorithm 2,
     * eliminating the need to enumerate O(n⁴) quartets as in the original wQFM.
     * 
     * @param book Current bookkeeping state containing taxa partitions and scoring data
     * @return true if the bipartition was improved, false if no improvement found
     */
    public static boolean oneInteration(BookKeepingPerLevelv2 book){
        
        // Cumulative gain tracking for determining optimal stopping point
        double cg = 0;                    // Current cumulative gain
        int maxCgIndex = -1;              // Index of maximum cumulative gain
        double maxCg = 0;                 // Maximum cumulative gain achieved

        // Special handling for singleton partitions (size 1)
        // This affects the stopping criteria as described in the paper
        boolean singletonPartition = book.taxas.getTaxonCountInPartition(0) == 1  || book.taxas.getTaxonCountInPartition(1) == 1;

        // Lock arrays to prevent taxa from being moved multiple times in one iteration
        // This is a core requirement of the FM algorithm to ensure convergence
        boolean[] rtLocked = new boolean[book.taxas.realTaxonCount];
        boolean[] dtLocked = new boolean[book.taxas.dummyTaxonCount];
        
        // Gain arrays for storing transfer benefits
        // rtGains[i][j] = gain from moving real taxon i to partition j
        // dtGains[i] = gain from moving dummy taxon i to opposite partition
        double[][] rtGains;
        double[] dtGains;

        // Track all swaps performed during this iteration for potential rollback
        ArrayList<Swap> swaps = new ArrayList<Swap>();

        /**
         * MAIN FM ITERATION LOOP
         * 
         * This loop implements the core FM strategy: repeatedly find and execute
         * the best taxon transfer until no valid transfers remain. The direct
         * scoring method from Algorithm 2 is used for gain calculations.
         */
        while(true){
            // Compute gains for all possible taxon transfers
            // This uses the direct scoring formulations from Section 2.5 of the paper
            rtGains = new double[book.taxas.realTaxonCount][2];
            dtGains = new double[book.taxas.dummyTaxonCount];
            
            // Calculate current bipartition score and all transfer gains
            // This implements the mathematical formulations from Section 2.5
            book.calculateScoreAndGains(rtGains, dtGains);

            // Find and execute the best taxon transfer
            var x = swapMax(book, rtGains, dtGains, rtLocked,dtLocked);
            
            if(x != null){
                // Record the swap for potential rollback
                swaps.add(x);
                
                double gain = x.gain;

                // Update cumulative gain
                cg += gain;

                // Track the best cumulative gain for optimal stopping
                // Special case for singleton partitions
                if(singletonPartition){
                    if(maxCgIndex == -1 ){
                        maxCg = cg;
                        maxCgIndex = swaps.size() - 1;
                    }
                }
                
                // Update maximum cumulative gain if current is better
                // The EPS threshold prevents numerical precision issues
                if(cg > maxCg && Math.abs(maxCg - cg) > EPS ){
                    maxCg = cg;
                    maxCgIndex = swaps.size() - 1;
                }

            }
            else{
                // No more valid swaps available, exit loop
                break;
            }
        }

        /**
         * OPTIMAL STOPPING POINT DETERMINATION
         * 
         * The FM algorithm may have performed swaps that decreased the overall
         * score after reaching the optimum. We now rollback to the point of
         * maximum cumulative gain to ensure we keep the best bipartition found.
         */
        if(maxCgIndex == -1){
            // No improvement found - rollback all swaps and indicate no progress
            if(swaps.size() != (book.taxas.realTaxonCount + book.taxas.dummyTaxonCount)){
                // Partial rollback: undo all swaps performed
                for(int i = swaps.size() - 1; i >= 0; --i){
                    var x = swaps.get(i);
                    if(x.isDummy){
                        book.taxas.swapPartitionDummyTaxon(x.index);
                    }
                    else{
                        book.taxas.swapPartitionRealTaxon(x.index);
                    }
                }
            }
            return false;
        }
        
        // Rollback swaps beyond the optimal point
        // Keep swaps up to maxCgIndex, undo the rest
        for(int i = swaps.size() - 1; i > maxCgIndex; --i){
            var x = swaps.get(i);
            book.swapTaxon(x.index, x.isDummy);
        }

        // Indicate that the bipartition was improved
        return true;
    }

}
