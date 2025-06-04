package src.ScoreCalculator;

import src.Config;
import src.Tree.Branch;

/**
 * NumSatCalculatorNodeEv2: Direct bipartition scoring from gene trees implementation.
 * 
 * This class implements Algorithm 2 from the paper: "Scoring a bipartition (used by FM 
 * algorithm to find the best bipartition)". It provides the mathematical formulations
 * described in Section 2.5 for computing bipartition scores directly from gene trees
 * without explicitly enumerating the O(n⁴) induced quartets.
 * 
 * Key mathematical concepts implemented:
 * 
 * 1. w(S^[g]) - Sum of weights of satisfied resolved quartets (Section 2.5.1)
 * 2. w(S^[g] ∪ V^[g] ∪ U^[g]) - Sum of weights of all relevant quartets (Section 2.5.2)
 * 3. w(U^[g]) - Sum of weights of unresolved quartets (Section 2.5.3)
 * 
 * The final score uses the restructured equation from the paper:
 * Score(A,B,G) = Σ_g (2w(S^[g]) - w(S^[g] ∪ V^[g] ∪ U^[g]) + w(U^[g]))
 * 
 * This reformulation enables efficient computation by avoiding direct calculation
 * of w(V^[g]) (violated quartets) while maintaining mathematical equivalence.
 */
public class NumSatCalculatorNodeEv2 implements NumSatCalculatorNode {
    
    // Gene tree branch structure representing components after removing internal nodes
    Branch[] branches;
    
    // Precomputed values for efficient quartet weight calculations
    double[] pairsBFromSingleBranch;    // w(PB_k^[g,u]) values for each branch
    double[][] sumPairsBranch;          // Sums of pair weights for each branch and partition
    double[] dummyTaxaWeightsIndividual; // Individual weights for dummy taxa
    double[] totalTaxa;                 // Total taxa counts per partition [A, B]
    double[][][] pairs;                 // Pair weights between branches [i][j][partition]
    double[] sumPairs;                  // Total pair weights per partition

    // Unresolved quartet handling
    double nonQuartets;                 // w(U^[g]) - sum of unresolved quartet weights
    double sumPairsBSingleBranch;       // Sum of single-branch B-partition pairs
    
    int nDummyTaxa;                     // Number of dummy taxa in current subproblem

    // Strategy pattern for different unresolved quartet calculation methods
    NonQuartCalculator nonQuartCalculator;

    /**
     * Interface for computing unresolved quartet contributions.
     * 
     * The paper describes two different approaches (Type A and Type B) for handling
     * unresolved quartets in polytomy nodes, as mentioned in Section 2.5.3.
     * This interface allows switching between these strategies.
     */
    public interface NonQuartCalculator{
        /**
         * Computes w(U^[g]) - the sum of weights of unresolved quartets.
         * @return Total weight of unresolved quartets in current gene tree
         */
        double calcNonQuartets();
        
        /**
         * Computes the change in unresolved quartet weights when transferring
         * taxa from one partition to another (used for gain calculations).
         * @param branchIndex Index of the branch affected by taxon transfer
         * @return Change in unresolved quartet weight
         */
        double changeAmount(int branchIndex);
    }

    /**
     * Type B unresolved quartet calculator.
     * 
     * This implements one variant of the unresolved quartet weight calculation
     * as described in Section 2.5.3 of the paper. The specific mathematical
     * formulation differs between Type A and Type B approaches.
     */
    // QUES: But, how is this correct?
    class NonQuartCalculatorB implements NonQuartCalculator{
        @Override
        public double calcNonQuartets(){
            double q = 0;
            // Calculate unresolved quartet weights using Type B formulation
            for(int i = 0; i < branches.length; ++i){
                for(int j = i + 1; j < branches.length; ++j){
                    // QUES: but how is this correct?
                    q += pairs[i][j][0] * (sumPairs[1] - pairs[i][j][1]);
                }
            }
            return q;
        }

        @Override
        public double changeAmount(int branchIndex){
            double q = 0;
            // Calculate change in unresolved quartet weights for gain computation
            for(int i = 0; i < branches.length; ++i){
                if(i == branchIndex) continue;
                int mni = branchIndex > i ? i : branchIndex;
                int mxi = branchIndex > i ? branchIndex : i;
                q += pairs[mni][mxi][0] * (sumPairs[1] - sumPairsBranch[branchIndex][1]);
                q += pairs[mni][mxi][1] * (sumPairs[0] - pairs[mni][mxi][0]);
            }
            return q;
        }
    }

    /**
     * Type A unresolved quartet calculator.
     * 
     * This implements the alternative variant of unresolved quartet weight calculation.
     * The mathematical formulation differs from Type B in how it handles the
     * interaction between branches in polytomy nodes.
     */
    class NonQuartCalculatorA implements NonQuartCalculator{
        @Override
        public double calcNonQuartets(){
            double q = 0;
            // Calculate unresolved quartet weights using Type A formulation
            for(int i = 0; i < branches.length; ++i){
                for(int j = i + 1; j < branches.length; ++j){
                    // this is just a restructuring, 
                    // similar to the supplementary section of original
                    q += pairs[i][j][0] * (sumPairs[1] - sumPairsBranch[i][1] - sumPairsBranch[j][1] + pairs[i][j][1]);
                }
            }
            return q;
        }

        @Override
        public double changeAmount(int branchIndex){
            double q = 0;
            // Calculate change in unresolved quartet weights for gain computation
            for(int i = 0; i < branches.length; ++i){
                if(i == branchIndex) continue;
                int mni = branchIndex > i ? i : branchIndex;
                int mxi = branchIndex > i ? branchIndex : i;
                // QUES: why this?
                // we are just moving one reaL raxon, right?
                // so, maybe only totalTaxaCount[branchIndex][0] etc...?
                q += pairs[mni][mxi][0] * (sumPairs[1] - sumPairsBranch[mni][1] - sumPairsBranch[mxi][1] + pairs[mni][mxi][1]);
                q += pairs[mni][mxi][1] * (sumPairs[0] - sumPairsBranch[mni][0] - sumPairsBranch[mxi][0] + pairs[mni][mxi][0]);
                
            }
            return q;
        }
    }

    // Dummy taxa partition assignment for weight calculations
    int[] dummyTaxaPartition;

    /**
     * Constructor initializes the scoring calculator for a specific gene tree node.
     * 
     * @param b Array of branches representing components after removing an internal node
     *          This corresponds to the C^[g,u] sets described in Section 2.5.1
     */
    public NumSatCalculatorNodeEv2(Branch[] b) {

        this.branches = b;
        this.totalTaxa = new double[2];

        // Initialize dummy taxa count from branch structure
        this.nDummyTaxa = b[0].dummyTaxaWeightsIndividual.length;

        // Initialize arrays for efficient quartet weight computation
        pairsBFromSingleBranch = new double[b.length];
        sumPairsBranch = new double[b.length][2];
        this.sumPairs = new double[2];
        this.sumPairsBSingleBranch = 0;
        this.pairs = new double[b.length][b.length][2];
        this.nonQuartets = 0;

        // Select unresolved quartet calculation strategy based on configuration
        if(Config.NON_QUARTET_TYPE == Config.NonQuartetType.A){
            this.nonQuartCalculator = new NonQuartCalculatorA();
        }
        else{
            this.nonQuartCalculator = new NonQuartCalculatorB();
        }
        
        // Initial calculation of unresolved quartet weights
        this.nonQuartets = this.nonQuartCalculator.calcNonQuartets();
    }

    /**
     * Initializes bookkeeping structures for bipartition scoring calculations.
     * 
     * This method implements the preprocessing required for efficient computation
     * of the mathematical formulations from Section 2.5. It precomputes various
     * pair weights and sums to enable O(1) access during score calculations.
     * 
     * Key computations performed:
     * 1. Initialize w(PA_{i,j}^[g,u]) and w(PB_{i,j}^[g,u]) calculations (Section 2.1.1)
     * 2. Precompute w(PB_k^[g,u]) for single-branch pairs (Section 2.1.2)
     * 3. Setup dummy taxa weight handling with normalization (Section 2.3)
     * 
     * @param dummyTaxaToPartitionMap Assignment of dummy taxa to partitions
     * @param totalTaxaA Total weight of taxa in partition A
     * @param totalTaxaB Total weight of taxa in partition B  
     * @param dummyTaxaWeightsIndividual Individual weights for dummy taxa
     * @param nDummyTaxa Number of dummy taxa in current subproblem
     */
    public void initBookkeeping(
        int[] dummyTaxaToPartitionMap,
        double totalTaxaA, 
        double totalTaxaB, 
        double[] dummyTaxaWeightsIndividual,
        int nDummyTaxa
    ){
        this.dummyTaxaPartition = dummyTaxaToPartitionMap;
        this.dummyTaxaWeightsIndividual = dummyTaxaWeightsIndividual;
        this.totalTaxa[0] = totalTaxaA;
        this.totalTaxa[1] = totalTaxaB;
        
        this.nDummyTaxa = nDummyTaxa;
        
        // Reset all precomputed values
        for(int i = 0; i < this.branches.length; ++i){
            this.pairsBFromSingleBranch[i] = 0;
            this.sumPairsBranch[i][0] = 0;
            this.sumPairsBranch[i][1] = 0;
        }
        this.sumPairs[0] = 0;
        this.sumPairs[1] = 0;
        this.sumPairsBSingleBranch = 0;

        var b = this.branches;
        
        /**
         * PRECOMPUTE INTER-BRANCH PAIR WEIGHTS
         * 
         * This implements the pair weight calculations described in Section 2.1.1.
         * For each pair of branches (i,j), compute:
         * - pairs[i][j][0] = w(PA_{i,j}^[g,u]) for partition A
         * - pairs[i][j][1] = w(PB_{i,j}^[g,u]) for partition B
         * 
         * The dummy taxa weight adjustments implement the constraint that
         * "no two taxa under the same dummy taxon" can participate in a quartet.
         */
        for(int i = 0; i < b.length; ++i){
            for(int j = i + 1; j < b.length; ++j){
                // Initial pair weights: product of branch total counts
                // pairs[i][j][0] = w(PA_{i,j}^[g,u]) for partition A
                // pairs[i][j][1] = w(PB_{i,j}^[g,u]) for partition B
                // here, since no weighting, just multiplying suffices
                this.pairs[i][j][0] = b[i].totalTaxaCounts[0] * b[j].totalTaxaCounts[0];
                this.pairs[i][j][1] = b[i].totalTaxaCounts[1] * b[j].totalTaxaCounts[1];
                
                // Subtract dummy taxa contributions to prevent two taxa under same dummy taxon
                // This implements the exclusion constraint from Section 2.5
                // all the dummy taxa are in one array kindof
                // for each dummy taxon, we find the partition it belongs to
                // and subtract from the pair product corresponding to that partition
                for(int k = 0; k < this.nDummyTaxa; ++k){
                    int partition = this.dummyTaxaPartition[k];
                    this.pairs[i][j][partition] -= b[i].dummyTaxaWeightsIndividual[k] * b[j].dummyTaxaWeightsIndividual[k];
                }
                
                // Accumulate sums for efficient access during scoring
                // here we store some further sums, 
                // one is, for each branch, sum over pair with others
                // other is, sum over all i,j
                sumPairsBranch[i][0] += this.pairs[i][j][0];
                sumPairsBranch[j][0] += this.pairs[i][j][0];
                sumPairsBranch[i][1] += this.pairs[i][j][1];
                sumPairsBranch[j][1] += this.pairs[i][j][1];

                this.sumPairs[0] += this.pairs[i][j][0];
                this.sumPairs[1] += this.pairs[i][j][1];
            }

            /**
             * PRECOMPUTE SINGLE-BRANCH PAIR WEIGHTS
             * 
             * This computes w(PB_k^[g,u]) as described in Section 2.1.2.
             * These represent pairs of taxa within the same branch that both
             * belong to partition B, implementing the formula:
             * 
             * w(PB_k^[g,u]) = (1/2) * (w(F_B^[g,u,k])² - w(R_B^[g,u,k]) - Σ_X w(X_R^[g,u,k])²)
             */

             // this is the first term
            pairsBFromSingleBranch[i] = b[i].totalTaxaCounts[1] * b[i].totalTaxaCounts[1];
            
            // Subtract dummy taxa squared weights (prevents {a,a} pairs)
            // this is the third term, subtracting for each dummy taxon
            for(int k = 0; k < this.nDummyTaxa; ++k){
                int partition = this.dummyTaxaPartition[k];
                if(partition == 1){
                    pairsBFromSingleBranch[i] -= b[i].dummyTaxaWeightsIndividual[k] * b[i].dummyTaxaWeightsIndividual[k];
                }
            }
            
            // Subtract real taxa weights (unit weights, so just count)
            // second term
            pairsBFromSingleBranch[i] -= b[i].realTaxaCounts[1];
            
            // Divide by 2 since we're counting unordered pairs
            pairsBFromSingleBranch[i] /= 2;
            this.sumPairsBSingleBranch += pairsBFromSingleBranch[i];
        }

        // Reinitialize unresolved quartet calculator with current configuration
        if(Config.NON_QUARTET_TYPE == Config.NonQuartetType.A){
            this.nonQuartCalculator = new NonQuartCalculatorA();
        }
        else{
            this.nonQuartCalculator = new NonQuartCalculatorB();
        }
        
        // Compute initial unresolved quartet weights
        this.nonQuartets = this.nonQuartCalculator.calcNonQuartets();
    }

    /**
     * Computes the bipartition score using the restructured equation from Section 2.5.
     * 
     * This implements the core scoring formula:
     * Score = Σ_g (2w(S^[g]) - w(S^[g] ∪ V^[g] ∪ U^[g]) + w(U^[g]))
     * 
     * The method computes only the contribution from this specific gene tree node,
     * which represents one internal node u in gene tree g. The total score is
     * obtained by summing over all internal nodes across all gene trees.
     * 
     * Key components:
     * 1. w(S^[g,u]) - satisfied quartets anchored at this node (computed via pairs)
     * 2. w(U^[g,u]) - unresolved quartets at this node (computed by nonQuartCalculator)
     * 
     * @return Contribution to total bipartition score from this gene tree node
     */
    @Override
    public double score(){
        double res = 0;
        
        /**
         * COMPUTE w(S^[g,u]) CONTRIBUTION
         * 
         * This implements the efficient formulation from Section 2.1.3:
         * w(S^[g,u]) = Σ_{i<j} w(PA_{i,j}^[g,u]) * (w(PB^[g,u]) - w(PB_i^[g,u]) - w(PB_j^[g,u]))
         * 
         * The formula is restructured for computational efficiency as:
         * -Σ_i w(PB_i^[g,u]) * Σ_j≠i w(PA_{i,j}^[g,u]) + w(PA^[g,u]) * w(PB^[g,u])
         */
        for(int i = 0; i < this.branches.length; ++i){
            res -=  pairsBFromSingleBranch[i] * sumPairsBranch[i][0];
        }
        res += this.sumPairs[0] * this.sumPairsBSingleBranch;

        /**
         * ADD UNRESOLVED QUARTET CONTRIBUTION
         * 
         * This adds (1/2) * w(U^[g,u]) to the score as described in Section 2.5.3.
         * The factor of 1/2 accounts for the way unresolved quartets are counted
         * in the overall scoring scheme.
         */
        res += (this.nonQuartets / 2);
        return res;
    }

    @Override
    public void swapRealTaxon(int branchIndex, int currPartition){
        
        this.nonQuartets -= this.nonQuartCalculator.changeAmount(branchIndex);
        for(int i = 0; i < this.branches.length; ++i){
            if(branchIndex == i){
                // calculating changes
                // for sumPairsBranch and sumPairs
                // note that, when, say a real taxon moves from A to B, it still is in the same branch
                // so, pair product weights where one branch is this branch, changes
                // as a whole, as the new one comes, the change is equal to the sum for all other branches
                // so we just do a subtraction
                this.sumPairsBranch[i][1 - currPartition] += (this.totalTaxa[1-currPartition] - this.branches[i].totalTaxaCounts[1-currPartition]);
                this.sumPairsBranch[i][currPartition] -= (this.totalTaxa[currPartition] - this.branches[i].totalTaxaCounts[currPartition]);

                this.sumPairs[1 - currPartition] += (this.totalTaxa[1-currPartition] - this.branches[i].totalTaxaCounts[1-currPartition]);
                this.sumPairs[currPartition] -= (this.totalTaxa[currPartition] - this.branches[i].totalTaxaCounts[currPartition]);

            }
            else{
                int mni = branchIndex > i ? i : branchIndex;
                int mxi = branchIndex > i ? branchIndex : i;

                // also, the pair weights for every pair like (branchIndex, i) will change too
                // like the following
                this.pairs[mni][mxi][currPartition] -= this.branches[i].totalTaxaCounts[currPartition];
                this.pairs[mni][mxi][1 - currPartition] += this.branches[i].totalTaxaCounts[1 - currPartition];

                // here, please note that, in the above case, we only handled the changes of the branchIndex
                // for which we are swapping the real taxon
                // but note that, sumPairsBranch will also change for other branches
                // How?
                // for some other branch i, consider all the pairs with it
                // among all of them, there is also the pair (i, branchIndex)
                // so, this pair will inc/dec and thus the following
                this.sumPairsBranch[i][currPartition] -= this.branches[i].totalTaxaCounts[currPartition];
                this.sumPairsBranch[i][1 - currPartition] += this.branches[i].totalTaxaCounts[1 - currPartition];

                // note that, here we no more change the sumPairs
                // because that would be redundant
                // the change that we did in the earlier case actually accounts for all as a whole
            }
            

        }

        if(currPartition == 1){
            // if we are moving from B to A
            // then, PB_k will decrease by how much?
            // that one taxon would pair up with the others in same branch common with parition 1(B)
            // so, taxa count - 1
            pairsBFromSingleBranch[branchIndex] -= this.branches[branchIndex].totalTaxaCounts[1] - 1;
            this.sumPairsBSingleBranch -= this.branches[branchIndex].totalTaxaCounts[1] - 1;                
        }
        else{
            // same as above
            // just where, we are moving from A to B
            pairsBFromSingleBranch[branchIndex] += this.branches[branchIndex].totalTaxaCounts[1];
            this.sumPairsBSingleBranch += this.branches[branchIndex].totalTaxaCounts[1];                
        }
        this.totalTaxa[currPartition] -= 1;
        this.totalTaxa[1 - currPartition] += 1;

        this.nonQuartets += this.nonQuartCalculator.changeAmount(branchIndex);

        branches[branchIndex].swapRealTaxa(currPartition);

    }
    @Override
    public void swapDummyTaxon(int dummyIndex, int currPartition){
        for(int i = 0; i < this.branches.length; ++i){
            double wi = this.branches[i].dummyTaxaWeightsIndividual[dummyIndex];

            for(int j = i + 1; j < this.branches.length; ++j){
                
                double wj = this.branches[j].dummyTaxaWeightsIndividual[dummyIndex];
                
                // here, see that we taking one from the dummy taxa for this branch
                // and another from the common taxa between the other partition and the other branch
                // note that, it already ensures that no two taxa are under same dummy taxon
                double inc = wi * this.branches[j].totalTaxaCounts[1 - currPartition] + wj * this.branches[i].totalTaxaCounts[1-currPartition];
                this.pairs[i][j][1 - currPartition] += inc;
                
                // while decreasing, 
                // we need to be a little more careful
                // now, inside the parentheses, wj is subtracted
                // to make sure that we do not include pairs under same this same dummy taxon
                double dec = wi * (this.branches[j].totalTaxaCounts[currPartition] - wj) + wj * (this.branches[i].totalTaxaCounts[currPartition] - wi);
                this.pairs[i][j][currPartition] -= dec;
                
                this.sumPairsBranch[i][1 - currPartition] += inc;
                this.sumPairsBranch[j][1 - currPartition] += inc;
                this.sumPairsBranch[i][currPartition] -= dec;
                this.sumPairsBranch[j][currPartition] -= dec;

                this.sumPairs[1 - currPartition] += inc;
                this.sumPairs[currPartition] -= dec;
            }

            if(currPartition == 1){
                // if moving from B to A
                // we are subtracting wi from inside to account for the condition
                // of not under same dummy taxon
                this.pairsBFromSingleBranch[i] -= (this.branches[i].totalTaxaCounts[1] - wi) * wi;
                this.sumPairsBSingleBranch -= (this.branches[i].totalTaxaCounts[1] - wi) * wi;
            }
            else{
                // if moving from A to B
                // in this case, no need to subtract
                // since, we are bringing this new dummy taxon onto here
                // so, no overlap is actually possible
                // since, these taxa under this dummy taxa did not exist there in the first place
                this.pairsBFromSingleBranch[i] += (this.branches[i].totalTaxaCounts[1]) * wi;
                this.sumPairsBSingleBranch += (this.branches[i].totalTaxaCounts[1]) * wi;                     
            }
        }

        this.totalTaxa[1 - currPartition] += this.dummyTaxaWeightsIndividual[dummyIndex];
        this.totalTaxa[currPartition] -= this.dummyTaxaWeightsIndividual[dummyIndex];


        for(int i = 0; i < this.branches.length; ++i){
            branches[i].swapDummyTaxon(dummyIndex, currPartition);
        }

        this.nonQuartets = this.nonQuartCalculator.calcNonQuartets();
        
    }

    // simulate gain of branches for both way
    // the idea here is to efficiently update the helping value arrays
    // and from there recalculate the score
    // but note that, the whole DP is not run again
    @Override
    public double[][] gainRealTaxa(double originalScore, double multiplier) {
        double[][] gainsOfBranches = new double[this.branches.length][2];
        for(int i = 0; i < branches.length; ++i){
            for(int p = 0; p < 2; ++p){
                if(this.branches[i].realTaxaCounts[p] > 0){
                    // we first swap
                    this.swapRealTaxon(i, p);
                    gainsOfBranches[i][p] = multiplier * (this.score() - originalScore);
                    // then swap back
                    this.swapRealTaxon(i, 1 - p);
                }
            }
        }
        return gainsOfBranches;
    }

    // similarly simulate gains of dummy taxa
    @Override
    public void gainDummyTaxa(double originalScore, double multiplier, double[] dummyTaxaGains) {
        for(int i = 0; i < this.nDummyTaxa; ++i){
            int currPartition = this.dummyTaxaPartition[i];
            this.swapDummyTaxon(i, currPartition);
            dummyTaxaGains[i] += multiplier * (this.score() - originalScore);
            this.swapDummyTaxon(i, 1 - currPartition);
        }
    }

    


    
}
