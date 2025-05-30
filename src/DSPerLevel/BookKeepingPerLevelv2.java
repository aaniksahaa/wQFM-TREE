package src.DSPerLevel;

import src.Config;
import src.Utility;
import src.InitialPartition.IMakePartition;
import src.PreProcessing.GeneTrees;
import src.Taxon.DummyTaxon;
import src.Taxon.RealTaxon;

/**
 * BookKeepingPerLevelv2: Central coordinator for Algorithm 2 scoring across all gene trees.
 * 
 * This class manages the data structures and orchestrates the scoring process described
 * in Algorithm 2 of the paper: "Scoring a bipartition (used by FM algorithm to find the
 * best bipartition)". It coordinates the application of the mathematical formulations
 * from Section 2.5 across all input gene trees.
 * 
 * Key responsibilities:
 * 1. Maintains BookKeepingPerTreev2 instances for each gene tree
 * 2. Aggregates scores from individual gene trees using the restructured equation:
 *    Score(A,B,G) = Σ_g (2w(S^[g]) - w(S^[g] ∪ V^[g] ∪ U^[g]) + w(U^[g]))
 * 3. Computes gain calculations for FM algorithm taxon transfers
 * 4. Handles dummy taxon weight normalization as described in Section 2.3
 * 5. Manages partition state changes during iterative bipartition refinement
 * 
 * This is a critical component that enables the O(n²k log n) time complexity per
 * FM iteration mentioned in Section 2.2, replacing the O(n⁴) complexity of original wQFM.
 */
public class BookKeepingPerLevelv2 {

    // Core data structures
    public final GeneTrees geneTrees;                    // Input gene trees
    public TaxaPerLevelWithPartition taxas;             // Current taxa partition state

    public boolean allowSingleton;

    // Array of per-tree scoring managers implementing Algorithm 2 for each gene tree
    // Each BookKeepingPerTreev2 handles the mathematical formulations from Section 2.5
    // for a single gene tree, computing w(S^[g]), w(S^[g] ∪ V^[g] ∪ U^[g]), and w(U^[g])
    public static BookKeepingPerTreev2[] bookKeepingPerTrees;

    /**
     * Helper method for taxon count with optional score normalization.
     * 
     * This handles the different weight assignment strategies mentioned in Section 2.3.
     * When normalization is disabled, flattened counts include all real taxa under
     * dummy taxa. When enabled, each dummy taxon contributes as a single unit.
     */
    public int getTotalTaxon(int p){
        if(Config.SCORE_NORMALIZATION_TYPE == Config.ScoreNormalizationType.NO_NORMALIZATION){
            return taxas.getTaxonCountFlattenedInPartition(p);
        }
        return taxas.getTaxonCountInPartition(p);
    }

    /**
     * Gets the individual weight for a dummy taxon based on normalization strategy.
     * 
     * This implements the weight normalization described in Section 2.3 of the paper.
     * With normalization disabled, dummy taxa weights equal their flattened real taxa count.
     * With normalization enabled, each dummy taxon gets unit weight (w(X) = 1).
     */
    public int getDummyTaxonIndiWeight(int index){
        if(Config.SCORE_NORMALIZATION_TYPE == Config.ScoreNormalizationType.NO_NORMALIZATION){
            return taxas.getFlattenedCount(index);
        }
        return 1;
    }

    /**
     * Constructor: Initialize bookkeeping structures for all gene trees.
     * 
     * Creates BookKeepingPerTreev2 instances for each gene tree, which will
     * implement the direct scoring formulations from Section 2.5.
     */
    public BookKeepingPerLevelv2(GeneTrees geneTrees){

        this.geneTrees = geneTrees;
        BookKeepingPerLevelv2.bookKeepingPerTrees = new BookKeepingPerTreev2[geneTrees.geneTrees.size()];

        initialBookKeeping();
    }

    /**
     * Resets bookkeeping structures for a new taxa partition configuration.
     * 
     * This is called when starting evaluation of a new bipartition candidate.
     * Each BookKeepingPerTreev2 reinitializes its data structures to compute
     * the scoring formulations from Algorithm 2 for the new partition.
     */
    public void resetBookKeeping(TaxaPerLevelWithPartition taxas){
        this.taxas = taxas;
        if(taxas.smallestUnit) return;
        for(int i = 0; i < geneTrees.geneTrees.size(); ++i){
            BookKeepingPerLevelv2.bookKeepingPerTrees[i].resetBookkeeping(taxas);
        }
    }


    /**
     * Initializes BookKeepingPerTreev2 instances for all gene trees.
     * 
     * Each instance will handle the direct scoring calculations from Algorithm 2
     * for its respective gene tree, implementing the mathematical formulations
     * from Sections 2.5.1, 2.5.2, and 2.5.3.
     */
    private void initialBookKeeping(){
        for(int i = 0; i < geneTrees.geneTrees.size(); ++i){
            BookKeepingPerLevelv2.bookKeepingPerTrees[i] = new BookKeepingPerTreev2(geneTrees.geneTrees.get(i));
        }
    }


    /**
     * Applies score normalization and computes final gains for FM algorithm.
     * 
     * This method handles the final step of gain calculation where the raw scores
     * computed from Algorithm 2 are normalized and converted to gain values for
     * the FM algorithm. The normalization ensures that the restructured scoring
     * equation produces mathematically equivalent results to the original formulation.
     */
    private double gainCalcFromSatWithNorm(double[][] realTaxaGains, double[] dummyTaxaGains, double totalScore){

        double[] dtTotals = new double[this.taxas.dummyTaxonCount];
        double currTotals = 0;
        
        // Aggregate total quartet counts from all gene trees for normalization
        for(var x : BookKeepingPerLevelv2.bookKeepingPerTrees){
            // totals[0] += x.totalQuartetsAfterSwap(1);
            // totals[1] += x.totalQuartetsAfterSwap(0);
            for(int i = 0;i < this.taxas.dummyTaxonCount; ++i){
                dtTotals[i] += x.totalQuartetsAfterDummySwap(i, 1 - taxas.inWhichPartitionDummyTaxonByIndex(i));
            }

            currTotals += x.totalQuartets();
        }

        // Compute normalized gains for real taxa transfers
        for(int i = 0; i < taxas.realTaxonCount; ++i){
            double total = 0;
            int partition = taxas.inWhichPartitionRealTaxonByIndex(i);

            for(var bookTree : BookKeepingPerLevelv2.bookKeepingPerTrees){
                total += bookTree.totalQuartetsAfterSwap(taxas.realTaxa[i].id, 1 - partition);
            }
            
            // System.out.println("total : " + total);

            // Utility.addArrayToFirst(realTaxaGains[i], this.gainsToAll);
            realTaxaGains[i][partition] += totalScore;
            realTaxaGains[i][partition] = Config.SCORE_EQN.scoreFromSatAndTotal(total, realTaxaGains[i][partition]);
        }

        // Compute normalized gains for dummy taxa transfers
        for(int i = 0; i < taxas.dummyTaxonCount; ++i){
            
            dummyTaxaGains[i] = Config.SCORE_EQN.scoreFromSatAndTotal(
                dtTotals[i],
                dummyTaxaGains[i] + totalScore
            );

        }
        
        // Compute normalized total score
        totalScore = Config.SCORE_EQN.scoreFromSatAndTotal(
            currTotals,
            totalScore
        );
        
        // Convert absolute scores to relative gains (change from current score)
        for (int i = 0; i < realTaxaGains.length; i++) {
            realTaxaGains[i][taxas.inWhichPartitionRealTaxonByIndex(i)] -= totalScore;
        }
        for (int i = 0; i < dummyTaxaGains.length; i++) {
            dummyTaxaGains[i] -= totalScore;
        }

        return totalScore;
    }

    
    /**
     * Computes the bipartition score using Algorithm 2 formulations.
     * 
     * This method implements the core scoring calculation from Algorithm 2:
     * Score(A,B,G) = Σ_g (2w(S^[g]) - w(S^[g] ∪ V^[g] ∪ U^[g]) + w(U^[g]))
     * 
     * For each gene tree g, the corresponding BookKeepingPerTreev2 computes:
     * - w(S^[g]): sum of weights of satisfied resolved quartets (Section 2.5.1)
     * - w(S^[g] ∪ V^[g] ∪ U^[g]): sum of weights of all relevant quartets (Section 2.5.2)
     * - w(U^[g]): sum of weights of unresolved quartets (Section 2.5.3)
     * 
     * The scores are aggregated across all gene trees and normalized according
     * to the configured scoring equation.
     */
    public double calculateScore(){
        double totalScore = 0;
        double currTotals = 0;

        // Aggregate scores from all gene trees using Algorithm 2
        for(var bookTree : BookKeepingPerLevelv2.bookKeepingPerTrees){
            for(var node : bookTree.nodesForScore){
                // System.out.println(node.index);
                // if(node.index == 11) {
                //     System.out.println("11 node");
                // }
                // System.out.println(node.index);
                
                // Compute score for this internal node using Algorithm 2 formulations
                double currScore = node.info.scoreCalculator.score();
                totalScore +=  currScore * node.frequency;
            }
            currTotals += bookTree.totalQuartets();
        }
        
        // Apply configured score normalization
        return Config.SCORE_EQN.scoreFromSatAndTotal(currTotals, totalScore);
    
    }
    

    /**
     * Computes bipartition score and gains for all possible taxon transfers.
     * 
     * This method implements the combined scoring and gain calculation required
     * by the FM algorithm. For each gene tree, it:
     * 1. Computes the current bipartition score using Algorithm 2
     * 2. Calculates gains for transferring each real taxon to the opposite partition
     * 3. Calculates gains for transferring each dummy taxon to the opposite partition
     * 4. Aggregates these values across the tree using post-order traversal
     * 
     * The gain calculations enable the FM algorithm to efficiently evaluate all
     * possible taxon transfers without recomputing the full score each time,
     * achieving the O(n²k log n) complexity mentioned in Section 2.2.
     */
    public double calculateScoreAndGains(double[][] realTaxaGains, double[] dummyTaxaGains){
        double totalScore = 0;
        
        // Process each gene tree separately and aggregate results
        for(var bookTree : BookKeepingPerLevelv2.bookKeepingPerTrees){
            double[] gainsToAll = new double[2];

            // Compute scores and gains at each internal node (Algorithm 2 application points)
            for(var node : bookTree.nodesForScore){
                // System.out.println("node index : " + node.index);
                
                // Apply Algorithm 2 scoring formulations at this node
                double score = node.info.scoreCalculator.score();
                
                // Compute gain arrays for real taxa transfers at this node
                var branchGains = node.info.scoreCalculator.gainRealTaxa(score, node.frequency);
                
                // Compute gains for dummy taxa transfers at this node
                node.info.scoreCalculator.gainDummyTaxa(score, node.frequency, dummyTaxaGains);
                
                score *= node.frequency;
                
                // System.out.println("Score at node: " + node.index + " " + score);
                // System.out.println(score);
                totalScore += score;
    
                // Distribute gains to child subtrees for bottom-up aggregation
                var childs = node.childs;
                for(int i = 0; i < childs.size(); ++i){
                    Utility.subArrayToFirst(branchGains[i], branchGains[childs.size()]);
                    childs.get(i).info.gainsForSubTree = branchGains[i];
                }
                Utility.addArrayToFirst(gainsToAll, branchGains[childs.size()]);
            }

            /**
             * Bottom-up gain aggregation through tree traversal.
             * 
             * This implements the efficient gain computation strategy where
             * gains computed at internal nodes are propagated down to leaf
             * nodes, enabling O(n) gain calculation per tree rather than O(n²).
             */
            for(int i = bookTree.nodesForGains.length - 1; i > -1; --i){
                var node = bookTree.nodesForGains[i];
                for (int j = 0; j < node.childs.size(); j++) {
                    var child = node.childs.get(j);
                    Utility.addArrayToFirst(child.info.gainsForSubTree, node.info.gainsForSubTree);
                }
                
                // Reset node gains after propagation
                node.info.gainsForSubTree[0] = 0;
                node.info.gainsForSubTree[1] = 0;
    
            }

            // Collect gains at leaf nodes (real taxa)
            for(var node : bookTree.geneTree.leaves){
                if(node == null) continue;
                if(taxas.isInRealTaxa(node.taxon.id)){
                    // if(!bookTree.geneTree.isTaxonPresent(node.taxon.id)){
                    //     System.out.println("Taxon not present");
                    //     System.exit(-1);
                    // }
                    Utility.addArrayToFirst(
                        realTaxaGains[taxas.getRealTaxonIndex(node.taxon.id)], 
                        node.info.gainsForSubTree
                    );
                }
                // Reset leaf gains after collection
                node.info.gainsForSubTree[0] = 0;
                node.info.gainsForSubTree[1] = 0;

            }
            
            // Add tree-wide gains to all real taxa present in this gene tree
            for(int i = 0; i < this.taxas.realTaxonCount; ++i){
                if(bookTree.geneTree.isTaxonPresent(this.taxas.realTaxa[i].id)){
                    Utility.addArrayToFirst(
                        realTaxaGains[i], 
                        gainsToAll
                    );
                }
            }

        }
        
        // Apply score normalization and finalize gain calculations
        totalScore = gainCalcFromSatWithNorm(realTaxaGains, dummyTaxaGains, totalScore);

        return totalScore;

    }

    
    /**
     * Executes a real taxon transfer between partitions.
     * 
     * This method updates the partition assignment for a real taxon and
     * propagates the change to all BookKeepingPerTreev2 instances so they
     * can update their Algorithm 2 data structures accordingly.
     */
    private void swapRealTaxon(int index){

        int partition = taxas.inWhichPartitionRealTaxonByIndex(index);
        taxas.swapPartitionRealTaxon(index);

        // System.out.println("swapping : " + taxas.realTaxa[index].label + " " + partition);
        
        // Update all per-tree data structures to reflect the taxon transfer
        for(var x : BookKeepingPerLevelv2.bookKeepingPerTrees){
            x.swapRealTaxon(taxas.realTaxa[index], partition);
        }
        
    }

    /**
     * Unified interface for taxon transfers (real or dummy).
     * 
     * This method is used by the FM algorithm to execute the best taxon
     * transfer found during each iteration of bipartition refinement.
     */
    public void swapTaxon(int index, boolean isDummy){
        if(isDummy) this.swapDummyTaxon(index);
        else this.swapRealTaxon(index);
    }

    /**
     * Executes a dummy taxon transfer between partitions.
     * 
     * This method updates the partition assignment for a dummy taxon and
     * propagates the change to all BookKeepingPerTreev2 instances. Dummy
     * taxon transfers affect the weight calculations described in Section 2.3.
     */
    private void swapDummyTaxon(int index){
        int partition = taxas.inWhichPartitionDummyTaxonByIndex(index);
        taxas.swapPartitionDummyTaxon(index);

        // Update all per-tree data structures to reflect the dummy taxon transfer
        for(var x : BookKeepingPerLevelv2.bookKeepingPerTrees){
            x.swapDummyTaxon(index, partition);
        }

    }

    /**
     * Implements the divide step of the divide-and-conquer algorithm.
     * 
     * This method creates two new subproblems from the current bipartition:
     * 1. Splits real and dummy taxa according to current partition assignments
     * 2. Creates new dummy taxa representing the opposite partitions
     * 3. Generates initial bipartitions for each subproblem using Algorithm 1
     * 4. Returns TaxaPerLevelWithPartition objects for recursive processing
     * 
     * This implements the divide step described in Section 2 of the paper,
     * including the dummy taxon tree structure from Section 2.3.
     */
    public TaxaPerLevelWithPartition[] divide(IMakePartition makePartition, boolean allowSingleton){
        RealTaxon[][] rts = new RealTaxon[2][];
        DummyTaxon[][] dts = new DummyTaxon[2][];

        // int[][] rtsPart = new int[2][];
        // int[][] dtsPart = new int[2][];

        // Allocate arrays for taxa in each partition
        for(int i = 0; i < 2; ++i){
            rts[i] = new RealTaxon[taxas.getRealTaxonCountInPartition(i)];
            dts[i] = new DummyTaxon[taxas.getDummyTaxonCountInPartition(i)];
            // var x = makePartition.makePartition(rts[i], dts[i]);
            // rtsPart[i] = x.realTaxonPartition;
            // dtsPart[i] = x.dummyTaxonPartition;
        }

        int[] index = new int[2];

        // Distribute real taxa to their assigned partitions
        for(var x : taxas.realTaxa){
            int part = taxas.inWhichPartition(x.id);
            rts[part][index[part]++] = x;
        }
        
        // Distribute dummy taxa to their assigned partitions
        index[0] = 0;
        index[1] = 0;
        int i = 0;
        for(var x : taxas.dummyTaxa){
            int part = taxas.inWhichPartitionDummyTaxonByIndex(i++);
            dts[part][index[part]++] = x;
        }

        /**
         * Create new dummy taxa for the divide step.
         * 
         * This implements the dummy taxon creation described in Section 2.3:
         * - X_1 represents taxa in partition B (added to partition A's subproblem)
         * - X_2 represents taxa in partition A (added to partition B's subproblem)
         * 
         * Each dummy taxon maintains a tree structure containing the real taxa
         * it represents, enabling proper weight normalization in subproblems.
         */
        DummyTaxon[] newDt = new DummyTaxon[2];
        
        // BookKeepingPerLevel[] bookKeepingPerLevels = new BookKeepingPerLevel[2];
        TaxaPerLevelWithPartition[] taxaPerLevelWithPartitions = new TaxaPerLevelWithPartition[2];
        for( i = 0; i < 2; ++i){
            // Create dummy taxon representing the opposite partition
            newDt[i] = new DummyTaxon(rts[1 - i], dts[1 - i]);
            
            // Add the new dummy taxon to the current partition's dummy taxa array
            DummyTaxon[] dtsWithNewDt = new DummyTaxon[dts[i].length + 1];
            for(int j = 0; j < dts[i].length; ++j){
                dtsWithNewDt[j] = dts[i][j];
            }
            dtsWithNewDt[dtsWithNewDt.length - 1] = newDt[i];

            // Generate initial bipartition for the subproblem using Algorithm 1
            if(rts[i].length + dtsWithNewDt.length > 3){

                var y = makePartition.makePartition(rts[i], dtsWithNewDt, true);
                taxaPerLevelWithPartitions[i] = new TaxaPerLevelWithPartition(
                    rts[i], dtsWithNewDt, 
                    y.realTaxonPartition, 
                    y.dummyTaxonPartition, 
                    this.geneTrees.realTaxaCount
                );
            }
            else{
                // Small subproblems get trivial partitions (will create star trees)
                taxaPerLevelWithPartitions[i] = new TaxaPerLevelWithPartition(
                    rts[i], dtsWithNewDt, 
                    null, null,
                    this.geneTrees.realTaxaCount
                );
            }
            
            // bookKeepingPerLevels[i] = new BookKeepingPerLevel(this.geneTrees,x);
        }

        return taxaPerLevelWithPartitions;
    }


}
