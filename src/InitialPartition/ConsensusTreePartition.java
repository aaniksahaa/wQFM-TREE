package src.InitialPartition;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;

import src.Config;
import src.DSPerLevel.BookKeepingPerLevelv2;
import src.DSPerLevel.TaxaPerLevelWithPartition;
import src.PreProcessing.GeneTrees;
import src.Taxon.DummyTaxon;
import src.Taxon.RealTaxon;
import src.Tree.Branch;
import src.Tree.Info;
import src.Tree.Tree;
import src.Tree.TreeNode;

/**
 * ConsensusTreePartition: Implementation of Algorithm 1 from the wQFM-TREE paper.
 * 
 * This class implements "Creation of initial bipartition for a taxa set which will be
 * further refined by the FM algorithm" as described in the supplementary material.
 * It corresponds to the consensus tree-based method mentioned in Section 2.4 of the
 * main paper for computing initial bipartitions.
 * 
 * The algorithm works by:
 * 1. Constructing a greedy consensus tree (majority rule extended tree) from gene trees
 * 2. Evaluating each edge in the consensus tree as a potential bipartition
 * 3. Handling dummy taxa placement using weighted assignment based on their real taxa
 * 4. Selecting the bipartition with the highest score using Algorithm 2 scoring
 * 
 * This approach eliminates the time-consuming step of sorting O(n⁴) weighted quartets
 * that was required in the original wQFM algorithm, as mentioned in Section 2.2.
 */
public class ConsensusTreePartition implements IMakePartition {

    // Core consensus tree structure used for initial bipartition generation
    Tree consTree;

    // Fallback random partitioner for edge cases
    RandPartition randPartition;

    // Algorithm state variables
    int taxonCount;
    GeneTrees trees;
    BookKeepingPerLevelv2 book;     // Scoring mechanism using Algorithm 2
    double score;

    // Distance matrix for consensus tree construction
    // Used to build the majority rule extended tree from gene trees
    public double[][] dist;
    Map<String, RealTaxon> taxaMap;
    

    /**
     * Constructor: Initialize consensus tree-based partition maker.
     * 
     * Reads the consensus tree from file and prepares the distance matrix.
     * The consensus tree is created using PAUP v4.0b10 as mentioned in Section 2.4
     * of the paper ("greedy consensus tree of the input gene trees, also known as
     * the majority rule extended tree").
     * 
     * @param filePath Path to consensus tree file in Newick format
     * @param taxaMap Mapping from taxon names to RealTaxon objects
     * @param trees Input gene trees for distance calculations
     */
    public ConsensusTreePartition(String filePath, Map<String, RealTaxon> taxaMap, GeneTrees trees) throws FileNotFoundException{
        Scanner scanner = new Scanner(new File(filePath));
        String line = scanner.nextLine();
        this.consTree = new Tree(line, taxaMap);
        this.taxonCount = taxaMap.size();
        scanner.close();

        randPartition = new RandPartition();
        this.trees = trees;
        dist = new double[this.taxonCount][this.taxonCount];
        this.taxaMap = taxaMap;

        // Build distance matrix from consensus tree for gene tree processing
        calculateDistanceMatrix();
        // for(var x : this.taxaMap.entrySet()){
        //     System.out.println(x.getKey() + " " + x.getValue().id);
        // }
        // printDistanceMatrix();
        
    }

    public void printDistanceMatrix() {
        
        for (int i = 0; i < taxonCount; i++) {
            for (int j = 0; j < taxonCount; j++) {
                System.out.print(dist[i][j] + " ");
            }
            System.out.println();
        }
    }

    

    /**
     * Recursively calculates distances from a node to all leaf nodes.
     * 
     * This is part of the consensus tree distance matrix computation used
     * for preprocessing gene trees, enabling efficient bipartition evaluation.
     */
    private void distToLeaves(TreeNode node, int from, double dist){
        if(node.isLeaf()){
            this.dist[from][node.taxon.id] = dist;
        }
        else{
            for(var child : node.childs){
                distToLeaves(child, from, dist + 1);
            }
        }
    }

    
    /**
     * Calculates distances from a specific taxon to all other taxa in consensus tree.
     * 
     * Used to build the distance matrix that will be used during gene tree
     * preprocessing and consensus-based initial bipartition generation.
     */
    private void distFromTaxon(RealTaxon rt){
        var leafNode = this.consTree.leaves[rt.id];
        var parent = leafNode.parent;
        int from = rt.id;
        double dist = 1;
        while(parent != null){
            for(var child : parent.childs){
                if(child != leafNode){
                    distToLeaves(child, from, dist + 1);
                }
            }
            leafNode = parent;
            parent = parent.parent;
            dist++;
        }
    }

    /**
     * Builds complete distance matrix from consensus tree.
     * 
     * This matrix is used throughout the algorithm for efficient distance
     * computations during gene tree processing and bipartition evaluation.
     */
    private void calculateDistanceMatrix(){
        for(var x : this.taxaMap.entrySet()){
            distFromTaxon(x.getValue());    
        }
    }



    /**
     * Assigns all taxa in a subtree to partition 1 (otherwise partition 0).
     * 
     * This implements the subtree-based bipartition assignment described in
     * Algorithm 1. Each edge in the consensus tree defines a potential bipartition
     * where one side of the edge goes to partition A and the other to partition B.
     */
    private void assignSubTreeToPartition(TreeNode node, int[] rtsp, Map<Integer, Integer> idToIndex){
        if(node.isLeaf()){
            if(idToIndex.containsKey(node.taxon.id)){
                rtsp[idToIndex.get(node.taxon.id)] = 1;
            }
        }
        else{
            for(var child : node.childs){
                assignSubTreeToPartition(child, rtsp,idToIndex);
            }
        }
    }


    /**
     * Evaluates the score of a bipartition defined by a consensus tree node.
     * 
     * This implements the candidate bipartition scoring described in Algorithm 1.
     * For each edge in the consensus tree, we:
     * 1. Create a bipartition (S_A, S_B) based on the edge
     * 2. Assign dummy taxa using weighted sums as described in Section 2.4
     * 3. Score the bipartition using Algorithm 2 (direct scoring from gene trees)
     * 
     * The method uses Algorithm 2's scoring mechanism to evaluate bipartitions
     * without enumerating quartets, implementing the formulations from Section 2.5.
     */
    double scoreForPartitionByNode(TreeNode node, RealTaxon[] rts, DummyTaxon[] dts, boolean allowSingleton){

        // Initialize partition assignments
        int[] rtsP = new int[rts.length];    // Real taxa partition assignments
        int[] dtsp = new int[dts.length];    // Dummy taxa partition assignments
        

        Map<Integer, Integer> idToIndex = new HashMap<>();
        int i = 0;
        for(var x : rts){
            idToIndex.put(x.id, i++);
        }
    
        // Assign real taxa based on consensus tree node subtree
        assignSubTreeToPartition(node, rtsP, idToIndex);

        // Assign dummy taxa based on weighted sums of their real taxa
        // This implements the weighted assignment strategy from Section 2.4:
        // "assign X to partition S_A if sum of weights of real taxa in X_R 
        // that belong to A is greater than sum that belong to B"
        for(i = 0; i < dts.length; ++i){
            if(node.info.branches[0].dummyTaxaWeightsIndividual[i] >= .5){
                dtsp[i] = 1;
            }
        }

        // Initialize scoring bookkeeping if needed
        if(this.book == null){
            TaxaPerLevelWithPartition taxas = new TaxaPerLevelWithPartition(rts, dts, rtsP, dtsp, this.taxonCount);
            this.book = new BookKeepingPerLevelv2(trees);
            this.book.resetBookKeeping(taxas);
        }
        else{
            // Optimization: only reset bookkeeping if partition changed significantly
            int rtCount = 0;
            int dtCount = 0;
            boolean changed = false;
            for(i = 0; i < rts.length; ++i){
                if(rtsP[i] != this.book.taxas.inWhichPartitionRealTaxonByIndex(i)){
                    changed = true;
                    rtCount++;
                    // book.swapTaxon(i, false);
                }
            }

            for(i = 0; i < dts.length; ++i){
                if(dtsp[i] != this.book.taxas.inWhichPartitionDummyTaxonByIndex(i)){
                    changed = true;
                    dtCount++;
                    // book.swapTaxon(i, true);
                }
            }
            if(!changed){
                // System.out.println("Not changed");
                return this.score;
            }
            else{
                // Heuristic: full reset vs incremental updates based on change magnitude
                if(rtCount + 9 * dtCount + 5 > dts.length){
                    TaxaPerLevelWithPartition taxas = new TaxaPerLevelWithPartition(rts, dts, rtsP, dtsp, this.taxonCount);
                    // this.book = new BookKeepingPerLevel(trees, taxas, allowSingleton);
                    this.book.resetBookKeeping(taxas);

                }
                else{
                    // Incremental updates through swapping
                    for(i = 0; i < rts.length; ++i){
                        if(rtsP[i] != this.book.taxas.inWhichPartitionRealTaxonByIndex(i)){
                            book.swapTaxon(i, false);
                        }
                    }

                    for(i = 0; i < dts.length; ++i){
                        if(dtsp[i] != this.book.taxas.inWhichPartitionDummyTaxonByIndex(i)){
                            book.swapTaxon(i, true);
                        }
                    }
                }
            }
            
        }
        
        // Calculate bipartition score using Algorithm 2 (direct scoring from gene trees)
        this.score = this.book.calculateScore();

        return this.score;
    }

    /**
     * Main method implementing Algorithm 1: Creation of initial bipartition.
     * 
     * This method implements the complete Algorithm 1 from the supplementary material:
     * 1. Construct greedy consensus tree from gene trees (done in constructor)
     * 2. For each edge e in consensus tree:
     *    a. Create candidate bipartition (A, B) defined by edge e
     *    b. Process real taxa: assign based on consensus tree partition
     *    c. Process dummy taxa: assign based on weighted sums of their real taxa
     *    d. Score candidate bipartition using Algorithm 2
     * 3. Return bipartition with highest score as initial bipartition
     * 
     * This eliminates the O(n⁴) quartet sorting step from original wQFM while
     * maintaining the quality of initial bipartitions for FM algorithm refinement.
     */
    @Override
    public MakePartitionReturnType makePartition(RealTaxon[] rts, DummyTaxon[] dts, boolean allowSingleton) {
        
        // Reset scoring bookkeeping for fresh evaluation
        this.book = null;

        // Weight assignment arrays for taxa in consensus tree
        double[] weight = new double[consTree.leavesCount];
        int[] inWhichDummyTaxa = new int[consTree.leavesCount];

        // Variables for tracking best bipartition candidate
        TreeNode minNode = null;
        double minDiff = 0;
        double maxScore = 0;

        boolean[] isRealTaxon = new boolean[consTree.leavesCount];

        // boolean allowSingleton = Config.ALLOW_SINGLETON;

        // Initialize weights for real taxa (unit weights as per Section 2.3)
        for(var x : rts){
            weight[x.id] = 1;
            isRealTaxon[x.id] = true;
        }

        int i = 0;

        /**
         * Weight normalization for dummy taxa as described in Section 2.3.
         * 
         * This implements the weight assignment w(a) from equation (1) in the paper:
         * - Real taxa not under dummy taxa get unit weight (w(a) = 1)
         * - Real taxa under dummy taxon X get fractional weights such that
         *   the total weight of X equals 1: w(X) = Σ_{a ∈ X_R} w(a) = 1
         * 
         * Two approaches are supported:
         * 1. NESTED: Non-uniform weights based on dummy taxon tree structure
         * 2. Default: Uniform weights (each real taxon under X gets weight 1/|X_R|)
         */
        for(var x : dts){
            if(Config.CONSENSUS_WEIGHT_TYPE == Config.ConsensusWeightType.NESTED){
                // Non-uniform normalization using tree structure (Section 2.3)
                x.calcDivCoeffs(Config.ScoreNormalizationType.NESTED_NORMALIZATION, weight, 1);
                for(var y : x.flattenedRealTaxa){
                    weight[y.id] = 1 / weight[y.id];
                }
            }
            else{
                // Uniform normalization: each real taxon gets equal weight
                double sz = x.flattenedTaxonCount;
                for(var y : x.flattenedRealTaxa){
                    weight[y.id] += 1. / sz;
                }
            }
            
            // Track which dummy taxon each real taxon belongs to
            for(var y : x.flattenedRealTaxa){
                inWhichDummyTaxa[y.id] = i;
            }
            ++i;
        }

        // if(Config.CONSENSUS_WEIGHT_TYPE == Config.ConsensusWeightType.NESTED){
        //     // for(i = 0; i < dts.length; ++i){
        //     //     if(weight[i] > 1) weight[i] = 1. / weight[i];
        //     // }
        // }        

        /**
         * Traverse consensus tree to evaluate all potential bipartitions.
         * 
         * For each internal node, we compute the weight distributions and
         * evaluate whether the bipartition defined by cutting the edge
         * above this node would be a valid candidate.
         */
        for(var node : this.consTree.topSortedNodes){
            // Initialize branch information for consensus tree node
            node.info = new Info();
            node.info.branches = new Branch[1];
            node.info.branches[0] = new Branch(dts.length, dts.length);

            var branch = node.info.branches[0];

            if(node.isLeaf()){
                // Leaf node: initialize with single taxon weight
                double w = weight[node.taxon.id];
                if(isRealTaxon[node.taxon.id]){
                    branch.realTaxaCounts[0] = 1;
                    branch.totalTaxaCounts[0] = 1;
                }
                else if(w != 0) {
                    branch.totalTaxaCounts[0] = w;
                    branch.dummyTaxaWeightsIndividual[inWhichDummyTaxa[node.taxon.id]] = w;
                }
            }
            else{
                // Internal node: aggregate weights from children
                for(var child : node.childs){

                    // Count taxa on each side of potential bipartition
                    int partASize = child.info.branches[0].realTaxaCounts[0];
                    int partBSize = rts.length - partASize;

                    // Aggregate dummy taxon weights and determine their assignments
                    for(int j = 0; j < dts.length; ++j){
                        branch.dummyTaxaWeightsIndividual[j] += child.info.branches[0].dummyTaxaWeightsIndividual[j];
                        
                        // Dummy taxon assignment based on weighted majority
                        // This implements the weighted assignment from Section 2.4
                        if(child.info.branches[0].dummyTaxaWeightsIndividual[j] >= .5){
                            partASize++;
                        }
                        else{
                            partBSize++;
                        }
                    }
                    branch.totalTaxaCounts[0] += child.info.branches[0].totalTaxaCounts[0];
                    branch.realTaxaCounts[0] += child.info.branches[0].realTaxaCounts[0];

                    /**
                     * Validate bipartition candidate and evaluate its quality.
                     * 
                     * A bipartition is valid if both sides have sufficient taxa:
                     * - Without singletons: both sides need > 1 taxon
                     * - With singletons: both sides need ≥ 1 taxon
                     * 
                     * Quality evaluation uses either:
                     * 1. Direct scoring via Algorithm 2 (if USE_SCORING_IN_CONSENSUS is true)
                     * 2. Balance-based heuristic (minimize size difference)
                     */
                    if((partASize > 1 && partBSize > 1) || (allowSingleton && partASize >= 1 && partBSize >= 1) ){
                        if(Config.USE_SCORING_IN_CONSENSUS){
                            // Use Algorithm 2 scoring for candidate evaluation
                            double score = scoreForPartitionByNode(child, rts, dts, allowSingleton);
                            if( minNode == null || score > maxScore){
                                maxScore = score;
                                minNode = child;
                            }
                        }
                        else{
                            // Use balance-based heuristic (minimize partition size difference)
                            double diff = Math.abs(rts.length + dts.length - child.info.branches[0].totalTaxaCounts[0]);
                            if(minNode == null || diff < minDiff){
                                minNode = child;
                                minDiff = diff;
                            }
                        }
                        // double diff = Math.abs(rts.length + dts.length - child.info.branches[0].totalTaxaCounts[0]);
                        // if(minNode == null || diff < minDiff){
                        //     minNode = child;
                        //     minDiff = diff;
                        // }
                        // else if(diff < minDiff){
                        //     minNode = child;
                        //     minDiff = diff;
                        // }
                    }

                }
            }
        }
        
        // Fallback to random partition if no valid consensus-based partition found
        if(minNode == null){
            System.out.println("Min Node null");
            return randPartition.makePartition(rts, dts, allowSingleton);
            // System.exit(-1);
        }
        
        // Create final bipartition assignment arrays
        // System.out.println("partition");
        int[] rtsP = new int[rts.length];     // Real taxa partition assignments
        int[] dtsp = new int[dts.length];     // Dummy taxa partition assignments

        Map<Integer, Integer> idToIndex = new HashMap<>();
        i = 0;
        for(var x : rts){
            idToIndex.put(x.id, i++);
        }
    
        // Assign real taxa based on selected consensus tree node
        assignSubTreeToPartition(minNode, rtsP, idToIndex);

        // Assign dummy taxa based on weighted majority assignment
        for(i = 0; i < dts.length; ++i){
            if(minNode.info.branches[0].dummyTaxaWeightsIndividual[i] >= .5){
                dtsp[i] = 1;
            }
        }

        return new MakePartitionReturnType(rtsP, dtsp);

    }
    
}
