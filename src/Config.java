package src;

// import src.ScoreCalculator.FractionSat;
import src.ScoreCalculator.SatSubVio;
import src.ScoreCalculator.ScoreEqn;

/**
 * Config: Central configuration class for wQFM-TREE algorithm variants.
 * 
 * This class defines the various algorithmic options and parameter settings
 * that control the behavior of the wQFM-TREE implementation. The configurations
 * correspond to different approaches described in the paper and supplementary
 * material for handling:
 * 
 * 1. Weight normalization strategies (Section 2.3)
 * 2. Consensus tree construction options (Algorithm 1)
 * 3. Scoring equation variations (Section 2.5)
 * 4. FM algorithm iteration control (Section 2.2)
 * 5. Unresolved quartet handling strategies (Section 2.5.3)
 * 
 * These settings allow experimentation with different algorithmic variants
 * to optimize performance for specific datasets and computational constraints.
 */
public class Config {
    
    /**
     * Weight normalization strategies for dummy taxa (Section 2.3).
     * 
     * Controls how weights are assigned to real taxa under dummy taxa:
     * - FLAT_NORMALIZATION: Uniform weights (each real taxon gets 1/|X_R|)
     * - NESTED_NORMALIZATION: Non-uniform weights based on dummy taxon tree structure
     * - NO_NORMALIZATION: No weight adjustment (flattened counting)
     * 
     * NESTED_NORMALIZATION corresponds to the weight assignment w(a) described
     * in equation (1) of the paper, ensuring Σ_{a ∈ X_R} w(a) = 1 for each
     * dummy taxon X while respecting the hierarchical tree structure.
     */
    public enum ScoreNormalizationType{
        FLAT_NORMALIZATION,
        NESTED_NORMALIZATION,
        NO_NORMALIZATION
    }
    
    // Default to nested normalization as described in Section 2.3
    public static ScoreNormalizationType SCORE_NORMALIZATION_TYPE = ScoreNormalizationType.NESTED_NORMALIZATION;

    /**
     * Consensus tree weight assignment strategies (Algorithm 1).
     * 
     * Controls how weights are computed during consensus tree-based initial
     * bipartition generation:
     * - FLAT: Uniform weight distribution
     * - NESTED: Hierarchical weight distribution based on dummy taxon structure
     * 
     * NESTED corresponds to the weighted assignment strategy described in
     * Section 2.4 where dummy taxa are assigned to partitions based on
     * weighted sums of their constituent real taxa.
     */
    public enum ConsensusWeightType{
        FLAT,
        NESTED
    }

    // Default to nested weighting for consensus tree operations
    public static ConsensusWeightType CONSENSUS_WEIGHT_TYPE = ConsensusWeightType.NESTED;

    /**
     * Controls whether singleton partitions are allowed during bipartition.
     * 
     * When false, both sides of a bipartition must contain > 1 taxon.
     * When true, partitions with single taxa are permitted, which can
     * improve algorithm flexibility but may affect tree quality.
     */
    public static boolean ALLOW_SINGLETON = false;
    
    /**
     * Enables level-based singleton threshold adjustment.
     * 
     * When true, the singleton threshold is adjusted based on recursion depth
     * to become more permissive in deeper levels of the divide-and-conquer
     * algorithm, helping with convergence on difficult subproblems.
     */
    public static boolean USE_LEVEL_BASED_SINGLETON_THRESHOLD = true;
    
    /**
     * Controls whether Algorithm 2 scoring is used in consensus tree evaluation.
     * 
     * When true, consensus tree bipartition candidates are evaluated using
     * the full Algorithm 2 scoring mechanism from Section 2.5.
     * When false, a simpler balance-based heuristic is used (minimize partition
     * size difference), which is faster but potentially less accurate.
     */
    public static boolean USE_SCORING_IN_CONSENSUS = true;

    /**
     * Threshold for singleton partition acceptance.
     * 
     * Controls the balance between allowing singleton partitions and maintaining
     * tree quality. Lower values are more restrictive.
     */
    public static double SINGLETON_THRESHOLD = .5; 
    
    /**
     * Multiplier for adjusting singleton threshold at deeper recursion levels.
     * 
     * Used when USE_LEVEL_BASED_SINGLETON_THRESHOLD is true to make the
     * algorithm more permissive at deeper levels of recursion.
     */
    public static double MAX_LEVEL_MULTIPLIER = .5;

    /**
     * Scoring equation implementation.
     * 
     * Controls which mathematical formulation is used for bipartition scoring:
     * - SatSubVio: "Satisfied minus Violated" scoring (default)
     * - FractionSat: Fractional satisfaction scoring
     * 
     * SatSubVio implements the restructured equation from Section 2.5:
     * Score(A,B,G) = Σ_g (2w(S^[g]) - w(S^[g] ∪ V^[g] ∪ U^[g]) + w(U^[g]))
     */
    public static ScoreEqn SCORE_EQN = new SatSubVio();
    // public static ScoreEqn SCORE_EQN = new FractionSat();

    /**
     * Maximum number of FM algorithm iterations per bipartition refinement.
     * 
     * Controls the convergence-quality tradeoff in the FM algorithm.
     * Higher values allow more refinement but increase computation time.
     * The O(n²k log n) time complexity mentioned in Section 2.2 is per iteration.
     */
    public static int MAX_ITERATION = 5;
    // public static boolean USE_MAX_DEPTH_MULTIPLIER = true;

    /**
     * Controls whether polytomies in gene trees are resolved before processing.
     * 
     * When true, non-binary nodes in gene trees are resolved using distance-based
     * heuristics before applying Algorithm 2 scoring. This ensures binary tree
     * structure required for efficient quartet evaluation.
     */
    public static boolean RESOLVE_POLYTOMY = false;

    /**
     * Strategy for handling unresolved quartets (Section 2.5.3).
     * 
     * Controls how unresolved quartets contribute to the scoring calculation:
     * - Type A: One approach to unresolved quartet weight computation
     * - Type B: Alternative approach with different weight handling
     * 
     * These correspond to different strategies for computing w(U^[g]) in the
     * restructured scoring equation from Section 2.5.
     */
    public enum NonQuartetType{
        A,B
    }

    // Default to Type A unresolved quartet handling
    public static NonQuartetType NON_QUARTET_TYPE = NonQuartetType.A;
}
