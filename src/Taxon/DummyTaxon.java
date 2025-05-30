package src.Taxon;

import src.Config;

/**
 * DummyTaxon: Implementation of dummy taxa tree structure from Section 2.3.
 * 
 * This class represents dummy taxa as described in the wQFM-TREE paper.
 * Dummy taxa are created during the divide step of the algorithm to represent
 * taxa from the opposite partition in each subproblem. They maintain a tree
 * structure of the real taxa they contain, enabling proper weight normalization.
 * 
 * Key features:
 * 1. Hierarchical structure: Dummy taxa can contain real taxa and other dummy taxa
 * 2. Weight normalization: Implements the constraint Σ_{a ∈ X_R} w(a) = 1
 * 3. Flattened access: Provides efficient access to all contained real taxa
 * 4. Nested levels: Tracks the depth of dummy taxon nesting for algorithm control
 * 
 * This implementation enables the divide-and-conquer approach where subproblems
 * maintain information about taxa from opposite partitions without explicitly
 * including them in the scoring calculations.
 */
public class DummyTaxon {

    // Global ID counter for unique dummy taxon identification
    private static int idCounter = 0;

    // Direct children of this dummy taxon
    public RealTaxon[] realTaxa;        // Real taxa directly under this dummy taxon
    public DummyTaxon[] dummyTaxa;      // Other dummy taxa directly under this dummy taxon
    
    // Flattened view of all real taxa under this dummy taxon (recursive)
    public RealTaxon[] flattenedRealTaxa;   // All real taxa contained under this dummy taxon

    // Size tracking for efficient operations
    public int taxonCount;              // Direct children count (realTaxa.length + dummyTaxa.length)
    public int realTaxonCount;          // Number of direct real taxa children
    public int flattenedTaxonCount;     // Total number of real taxa under this dummy taxon (recursive)

    // Identification and structure
    public int id;                      // Unique identifier for this dummy taxon
    public int nestedLevel;             // Depth of nesting (how many dummy taxa levels deep)


    /**
     * Constructor: Creates a new dummy taxon from real and dummy taxa.
     * 
     * This constructor implements the dummy taxon creation described in Section 2.3.
     * It builds the tree structure that will be used for weight normalization
     * and maintains the constraint that each dummy taxon represents exactly
     * one unit of weight: w(X) = Σ_{a ∈ X_R} w(a) = 1.
     * 
     * @param rts Real taxa to be directly contained in this dummy taxon
     * @param dts Dummy taxa to be directly contained in this dummy taxon
     */
    public DummyTaxon(RealTaxon[] rts, DummyTaxon[] dts){
        this.realTaxonCount = rts.length;
        this.taxonCount = rts.length + dts.length;
        this.realTaxa = rts;
        this.dummyTaxa = dts;

        this.nestedLevel = 0;

        // Calculate flattened count and maximum nesting level from contained dummy taxa
        for(var x : dts){
            this.flattenedTaxonCount += x.flattenedTaxonCount;
            this.nestedLevel = Math.max(this.nestedLevel, x.nestedLevel);
        }
        
        // Build flattened array of all contained real taxa for efficient access
        this.flattenedTaxonCount += this.realTaxonCount;
        this.flattenedRealTaxa = new RealTaxon[this.flattenedTaxonCount];
        
        int i = 0;
        
        // Add direct real taxa first
        for(var x : realTaxa){
            this.flattenedRealTaxa[i] = x;
            i++;
        }
        
        // Add all real taxa from contained dummy taxa (flattened)
        for(var x : dummyTaxa){
            for(var y : x.flattenedRealTaxa){
                this.flattenedRealTaxa[i] = y;
                i++;
            }
        }

        // Assign unique ID and increment nesting level
        this.id = idCounter++;
        this.nestedLevel += 1;

    }

    /**
     * Calculates weight normalization coefficients as described in Section 2.3.
     * 
     * This method implements the weight assignment w(a) from equation (1) in the paper.
     * It ensures that for each dummy taxon X, the sum of weights of real taxa under X
     * equals 1: Σ_{a ∈ X_R} w(a) = 1, while supporting different normalization strategies.
     * 
     * The method recursively traverses the dummy taxon tree structure to assign
     * appropriate weight coefficients to all contained real taxa.
     * 
     * @param normalizationType Strategy for weight normalization (FLAT, NESTED, or NONE)
     * @param coeffs Array to store weight coefficients indexed by real taxon ID
     * @param multiplier Current weight multiplier from parent dummy taxa
     */
    public void calcDivCoeffs(Config.ScoreNormalizationType normalizationType, double[] coeffs, double multiplier){
        if(normalizationType == Config.ScoreNormalizationType.NO_NORMALIZATION){
            // No normalization: all real taxa get unit weight
            for(var x : this.flattenedRealTaxa)
                coeffs[x.id] = 1;
        }
        else if(normalizationType == Config.ScoreNormalizationType.FLAT_NORMALIZATION){
            // Flat normalization: uniform weight distribution
            // Each real taxon gets weight 1/|X_R| where |X_R| is total real taxa under this dummy taxon
            for(var x : this.flattenedRealTaxa)
                coeffs[x.id] = this.flattenedTaxonCount;
        }
        else if(normalizationType == Config.ScoreNormalizationType.NESTED_NORMALIZATION){
            /**
             * Nested normalization: Weight distribution respecting tree structure.
             * 
             * This implements the hierarchical weight assignment described in Section 2.3.
             * The weight is distributed proportionally among direct children, then
             * recursively distributed within each dummy taxon child.
             * 
             * For this dummy taxon with children C = {real taxa ∪ dummy taxa}:
             * - Each child gets weight multiplier / |C|
             * - Real taxa receive this weight directly
             * - Dummy taxa recursively distribute their weight among their children
             */
            double sz = this.realTaxonCount + this.dummyTaxa.length;
            
            // Assign weight coefficients to direct real taxa children
            for(var x : this.realTaxa)
                coeffs[x.id] = sz * multiplier;
            
            // Recursively assign weights within dummy taxa children
            for(var x : this.dummyTaxa){
                x.calcDivCoeffs(normalizationType, coeffs, sz * multiplier);
            }
        }
    }
}
