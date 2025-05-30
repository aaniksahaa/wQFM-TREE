package src.DSPerLevel;

import java.util.ArrayList;

import src.Config;
import src.Taxon.DummyTaxon;
import src.Taxon.RealTaxon;
import src.Tree.Info;
import src.Tree.Tree;
import src.Tree.TreeNode;

/**
 * TaxaPerLevelWithPartition: Manages taxa partition state during wQFM-TREE execution.
 * 
 * This class encapsulates the current bipartition state for a set of real and dummy taxa
 * during the divide-and-conquer algorithm execution. It handles the complex weight
 * normalization scheme described in Section 2.3 of the paper and provides efficient
 * access methods for Algorithm 2 scoring.
 * 
 * Key responsibilities:
 * 1. Maintains current partition assignments for real and dummy taxa
 * 2. Implements weight normalization as described in Section 2.3:
 *    - Real taxa not under dummy taxa: w(a) = 1
 *    - Real taxa under dummy taxon X: w(a) such that Σ_{a ∈ X_R} w(a) = 1
 * 3. Provides efficient lookup methods for Algorithm 2 scoring operations
 * 4. Handles taxon transfers for FM algorithm iteration
 * 5. Creates star trees for base cases (≤3 taxa)
 * 
 * The partition state represents a candidate bipartition (A, B) where each taxon
 * is assigned to either partition 0 or partition 1, corresponding to sides A and B
 * of the bipartition evaluation in Algorithm 2.
 */
public class TaxaPerLevelWithPartition {

    public final boolean LEFT = false;
    public final boolean RIGHT = true;
    

    // Total number of real taxa in the original problem (for array sizing)
    public final int allRealTaxaCount;
    public final int divCoeffsType = 0;
    
    // Core taxa arrays for this subproblem
    public RealTaxon[] realTaxa;        // Real taxa in current subproblem
    public DummyTaxon[] dummyTaxa;      // Dummy taxa in current subproblem
    
    // Partition assignment arrays (0 or 1 for each taxon)
    public int[] realTaxonPartition;    // realTaxa[i] assigned to partition realTaxonPartition[i]
    public int[] dummyTaxonPartition;   // dummyTaxa[i] assigned to partition dummyTaxonPartition[i]

    public int realTaxonCount;
    public int dummyTaxonCount;

    // Fast lookup arrays indexed by global real taxon IDs
    private boolean[] isInRealTaxa;     // isInRealTaxa[id] = true if real taxon id is in this subproblem
    private boolean[] isInDummyTaxa;    // isInDummyTaxa[id] = true if real taxon id is under some dummy taxon
    private double[] coeffs;            // Weight normalization coefficients from Section 2.3
    
    private int[] inWhichPartition;     // inWhichPartition[id] = partition of real taxon id
    private int[] inWhichDummyTaxa;     // inWhichDummyTaxa[id] = index of dummy taxon containing real taxon id
    private int[] realTaxonIndex;       // realTaxonIndex[id] = index in realTaxa array for real taxon id

    // Partition size tracking for efficient Algorithm 2 scoring
    private int[] taxonCountsInPartitions;              // Total taxa count per partition
    private int[] realTaxonCountsInPartitions;          // Real taxa count per partition
    private int[] dummyTaxonCountsInPartitions;         // Dummy taxa count per partition
    private int[] dummyTaxonCountsFlattenedInPartitions;// Flattened real taxa under dummy taxa per partition

    // Flag indicating if this is a base case (≤3 taxa total)
    public boolean smallestUnit;

    /**
     * Constructor: Initialize taxa partition state.
     * 
     * Creates a new partition state from the given taxa arrays and partition assignments.
     * Implements the weight normalization scheme from Section 2.3 of the paper and
     * sets up efficient lookup structures for Algorithm 2 scoring operations.
     * 
     * @param rts Real taxa in this subproblem
     * @param dts Dummy taxa in this subproblem
     * @param rtp Partition assignments for real taxa (0 or 1 for each)
     * @param dtp Partition assignments for dummy taxa (0 or 1 for each)
     * @param rtc Total number of real taxa in original problem (for array sizing)
     */
    public TaxaPerLevelWithPartition( RealTaxon[] rts, DummyTaxon[] dts, int[] rtp, int[] dtp, int rtc){
        this.realTaxa = rts;
        this.dummyTaxa = dts;
        this.realTaxonPartition = rtp;
        this.dummyTaxonPartition = dtp;
        this.realTaxonCount = rts.length;
        this.dummyTaxonCount = dts.length;
        this.allRealTaxaCount = rtc;

        // Check for base case: subproblems with ≤3 taxa create star trees directly
        if(rts.length + dts.length < 4){
            this.smallestUnit = true;
            return;
        }
        this.smallestUnit = false;
        
        // Initialize lookup arrays (sized for all possible real taxa)
        this.isInRealTaxa = new boolean[this.allRealTaxaCount];
        this.coeffs = new double[this.allRealTaxaCount];           // Weight normalization coefficients
        this.realTaxonIndex = new int[this.allRealTaxaCount];
        this.taxonCountsInPartitions = new int[2];
        this.inWhichDummyTaxa = new int[this.allRealTaxaCount];
        this.isInDummyTaxa = new boolean[this.allRealTaxaCount];
        this.realTaxonCountsInPartitions = new int[2];
        this.dummyTaxonCountsInPartitions = new int[2];

        this.dummyTaxonCountsFlattenedInPartitions = new int[2];

        int i = 0;

        /**
         * Process real taxa: Initialize lookup structures and count partitions.
         * 
         * Real taxa not under dummy taxa get unit weight (w(a) = 1) as described
         * in Section 2.3. The partition counts are used for Algorithm 2 scoring.
         */
        for(var x : realTaxa){
            isInRealTaxa[x.id] = true;
            coeffs[x.id] = 1.0;                                     // Unit weight for real taxa
            taxonCountsInPartitions[realTaxonPartition[i]]++;
            realTaxonIndex[x.id] = i++;
        }
        this.realTaxonCountsInPartitions[0] = taxonCountsInPartitions[0];
        this.realTaxonCountsInPartitions[1] = taxonCountsInPartitions[1];

        /**
         * Process dummy taxa: Apply weight normalization and track containment.
         * 
         * This implements the weight normalization described in Section 2.3:
         * For each dummy taxon X, the sum of weights of real taxa under X equals 1:
         * Σ_{a ∈ X_R} w(a) = 1
         * 
         * The calcDivCoeffs method computes the appropriate weight coefficients
         * based on the dummy taxon tree structure and configured normalization type.
         */
        i = 0;
        for(var x : dummyTaxa){
            // Mark all real taxa under this dummy taxon
            for(var y : x.flattenedRealTaxa){
                isInDummyTaxa[y.id] = true;
                inWhichDummyTaxa[y.id] = i;
            }
            
            // Apply weight normalization to real taxa under this dummy taxon
            x.calcDivCoeffs(Config.SCORE_NORMALIZATION_TYPE, coeffs, 1.);
            
            // Update partition counts
            taxonCountsInPartitions[dummyTaxonPartition[i]]++;
            this.dummyTaxonCountsInPartitions[dummyTaxonPartition[i]]++;
            this.dummyTaxonCountsFlattenedInPartitions[dummyTaxonPartition[i]] += x.flattenedTaxonCount;
            
            ++i;
        }
        
        // Build partition assignment lookup for all real taxa
        i = 0;
        this.inWhichPartition = new int[this.allRealTaxaCount];
        for(var x : realTaxa){
            this.inWhichPartition[x.id] = realTaxonPartition[i++];
        }
        
        // Assign partitions to real taxa under dummy taxa
        i = 0;
        for(var x : dummyTaxa){
            for(var y : x.flattenedRealTaxa){
                this.inWhichPartition[y.id] = dummyTaxonPartition[i];
            }
            ++i;
        }
    }

    // Fast lookup methods for Algorithm 2 scoring

    /**
     * Checks if a real taxon (by global ID) is part of this subproblem.
     */
    public boolean isInRealTaxa(int realTaxonId){
        return isInRealTaxa[realTaxonId];
    }

    /**
     * Gets the normalized weight for a real taxon.
     * 
     * Returns the weight w(a) as described in Section 2.3. For real taxa not
     * under dummy taxa, this is 1. For real taxa under dummy taxon X, this
     * is the fractional weight such that Σ_{a ∈ X_R} w(a) = 1.
     */
    public double getWeight(int realTaxonId){
        return 1. / this.coeffs[realTaxonId];
    }

    /**
     * Gets the partition assignment (0 or 1) for a real taxon.
     */
    public int inWhichPartition(int realTaxaId){
        return inWhichPartition[realTaxaId];
    }

    /**
     * Checks if a real taxon is under some dummy taxon in this subproblem.
     */
    public boolean isInDummyTaxa(int realTaxonId){
        return isInDummyTaxa[realTaxonId];
    }
    
    /**
     * Gets the index of the dummy taxon containing a real taxon.
     */
    public int inWhichDummyTaxa(int realTaxonId){
        return this.inWhichDummyTaxa[realTaxonId];
    }
    
    /**
     * Gets the local array index for a real taxon in the realTaxa array.
     */
    public int getRealTaxonIndex(int realTaxonId){
        return this.realTaxonIndex[realTaxonId];
    }
    
    /**
     * Gets the total number of taxa (real + dummy) in a partition.
     */
    public int getTaxonCountInPartition(int partition){
        return this.taxonCountsInPartitions[partition];
    }

    // Partition assignment accessors by local index

    /**
     * Gets partition assignment for a real taxon by its local array index.
     */
    public int inWhichPartitionRealTaxonByIndex(int index){
        return this.realTaxonPartition[index];
    }
    
    /**
     * Gets partition assignment for a dummy taxon by its local array index.
     */
    public int inWhichPartitionDummyTaxonByIndex(int index){
        return this.dummyTaxonPartition[index];
    }
    
    /**
     * Gets the number of real taxa in a partition.
     */
    public int getRealTaxonCountInPartition(int partition){
        return this.realTaxonCountsInPartitions[partition];
    }
    
    /**
     * Gets the number of dummy taxa in a partition.
     */
    public int getDummyTaxonCountInPartition(int partition){
        return this.dummyTaxonCountsInPartitions[partition];
    }

    /**
     * Gets the flattened count of real taxa under dummy taxa in a partition.
     * 
     * This counts all real taxa that are contained within dummy taxa assigned
     * to the given partition, used for score normalization calculations.
     */
    public int getDummyTaxonCountFlattenedInPartition(int partition){
        return this.dummyTaxonCountsFlattenedInPartitions[partition];
    }
    
    /**
     * Gets the total flattened taxon count in a partition.
     * 
     * This includes real taxa plus all real taxa under dummy taxa, providing
     * the total effective taxon count for scoring purposes.
     */
    public int getTaxonCountFlattenedInPartition(int partition){
        return (this.realTaxonCountsInPartitions[partition] + this.dummyTaxonCountsFlattenedInPartitions[partition]);
    }

    /**
     * Gets the flattened real taxa count for a dummy taxon.
     */
    public int getFlattenedCount(int index){
        return this.dummyTaxa[index].flattenedTaxonCount;
    }

    /**
     * Transfers a real taxon between partitions.
     * 
     * This method is used by the FM algorithm to execute taxon transfers during
     * bipartition refinement. It updates all relevant counts and lookup structures
     * to maintain consistency for Algorithm 2 scoring.
     */
    public void swapPartitionRealTaxon(int index){
        int currPartition = this.realTaxonPartition[index];
        int switchedPartition = (int) (1 - currPartition);

        // Update partition assignment
        this.realTaxonPartition[index] = switchedPartition;
        this.inWhichPartition[this.realTaxa[index].id] = switchedPartition;
        
        // Update partition counts
        this.taxonCountsInPartitions[currPartition]--;
        this.taxonCountsInPartitions[switchedPartition]++;
        this.realTaxonCountsInPartitions[currPartition]--;
        this.realTaxonCountsInPartitions[switchedPartition]++;
    }
    

    /**
     * Transfers a dummy taxon between partitions.
     * 
     * This method transfers a dummy taxon and all real taxa under it to the
     * opposite partition. It updates all relevant counts and assignments to
     * maintain consistency for Algorithm 2 scoring operations.
     */
    public void swapPartitionDummyTaxon(int index){
        int currPartition = this.dummyTaxonPartition[index];
        int switchedPartition = (1 - currPartition);

        // Update dummy taxon partition assignment
        this.dummyTaxonPartition[index] = switchedPartition;
        
        // Update partition assignments for all real taxa under this dummy taxon
        for(var x : this.dummyTaxa[index].flattenedRealTaxa){
            this.inWhichPartition[x.id] = switchedPartition;
        }
        
        // Update partition counts
        this.taxonCountsInPartitions[currPartition]--;
        this.taxonCountsInPartitions[switchedPartition]++;
        this.dummyTaxonCountsInPartitions[currPartition]--;
        this.dummyTaxonCountsInPartitions[switchedPartition]++;

        // Update flattened counts (real taxa under dummy taxa)
        this.dummyTaxonCountsFlattenedInPartitions[currPartition] -= this.dummyTaxa[index].flattenedTaxonCount;
        this.dummyTaxonCountsFlattenedInPartitions[switchedPartition] += this.dummyTaxa[index].flattenedTaxonCount;

    }


    /**
     * Creates a star tree for base cases (≤3 taxa).
     * 
     * This method handles the base case of the divide-and-conquer algorithm
     * when a subproblem has ≤3 taxa. It creates a star tree (single internal
     * node with all taxa as leaves) which represents the trivial solution
     * for small subproblems.
     * 
     * The star tree structure ensures that the conquer phase can properly
     * combine subproblem solutions using the tree grafting operations
     * described in the paper.
     */
    public Tree createStar(){
        if(!smallestUnit){
            System.out.println("Create Star should be called only on smallest unit\n");
            System.exit(-1);
        }

        Tree t = new Tree();
        ArrayList<TreeNode> childs = new ArrayList<>();
        
        // Add real taxa as leaf nodes
        for(var x : this.realTaxa){
            childs.add(t.addLeaf(x).setInfo(new Info(-1)));
        }

        // Add dummy taxa as leaf nodes (with special dummy taxon info)
        for(var x : this.dummyTaxa){
            childs.add(t.addLeaf(null).setInfo(new Info(x.id)));
        }

        // Create root internal node connecting all leaves (star structure)
        t.root = t.addInternalNode(childs).setInfo(new Info(-1));

        return t;
    }

}
