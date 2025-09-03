package src.PreProcessing;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;

import src.Config;
import src.Taxon.RealTaxon;
import src.Tree.Tree;
import src.Tree.TreeNode;

/**
 * GeneTrees: Preprocessing and management of input gene trees for wQFM-TREE.
 * 
 * This class handles the critical first step of the wQFM-TREE algorithm: preprocessing
 * the input gene trees for efficient use in Algorithm 2 scoring. It performs several
 * key operations described in the paper:
 * 
 * 1. Taxa name extraction and standardization across all gene trees
 * 2. Gene tree parsing from Newick format with consistent taxon mapping
 * 3. Optional polytomy resolution for binary tree structure requirement
 * 4. Tripartition frequency calculation for consensus tree construction
 * 5. Data structure initialization for efficient Algorithm 2 operations
 * 
 * The preprocessing ensures that all gene trees use consistent taxon IDs and
 * are structured appropriately for the O(n²k log n) scoring complexity described
 * in Section 2.2. The tripartition data supports Algorithm 1 for initial
 * bipartition generation.
 */
public class GeneTrees {

    // Core data structures
    public ArrayList<Tree> geneTrees;               // Parsed and preprocessed gene trees
    public String[] taxonIdToLabel;                 // ID to label mapping for output
    public RealTaxon[] taxa;                        // Array of all real taxa (indexed by ID)
    public Map<String, TreeNode> triPartitions;     // Tripartition frequency data for consensus
    public Map<String, RealTaxon> taxaMap;          // Label to RealTaxon mapping
    public int realTaxaCount;                       // Total number of real taxa
    public String path;                             // Input file path

    /**
     * Extracts taxon names from a single Newick tree string with support values and branch lengths.
     * 
     * Enhanced parser that handles both simple and extended Newick formats:
     * - Simple: ((1,2),3);
     * - Extended: ((A:1,B:2)0.98:3,E:2);
     * - OR LIKE (A:1,B:2,(F:1,(E:3,(C:2,D:5)98:4)88:1)87:3); where the support values are integers 
     */
    // private void parseTaxa(String newickLine, Set<String> taxaSet){
    //     newickLine = newickLine.replaceAll("\\s", "");
    //     int n = newickLine.length();
    //     int i = 0;
    
    //     while(i < n){
    //         char curr = newickLine.charAt(i);
    //         if(curr == '(' || curr == ')' || curr == ',' || curr == ';'){
    //             // Structural characters - skip
    //             i++;
    //         }
    //         else if(curr == ':'){
    //             // Branch length marker - skip all following numbers
    //             i++;
    //             while(i < n && (Character.isDigit(newickLine.charAt(i)) || newickLine.charAt(i) == '.')) {
    //                 i++;
    //             }
    //         }
    //         else{
    //             // Extract potential taxon name or support value
    //             StringBuilder token = new StringBuilder();
    //             int j = i;
                
    //             // Extract token until we hit a structural character or branch marker
    //             while(j < n){
    //                 char curr_j = newickLine.charAt(j);
    //                 if(curr_j == ':' || curr_j == ')' || curr_j == ',' || curr_j == ';' || curr_j == '('){
    //                     break;
    //                 }
    //                 token.append(curr_j);
    //                 j++;
    //             }
                
    //             String label = token.toString();
                
    //             // Add as taxon if it's non-empty and not a support value
    //             if(label.length() > 0 && !isSupportValue(label, j, newickLine)){
    //                 taxaSet.add(label);
    //             }
                
    //             // Move to the position after this token
    //             i = j;
    //         }
    //     }
    // }
    
    /**
     * Helper method to check if a string represents a support value.
     * 
     * Support values can be:
     * - Floating point between 0 and 1: "0.95", "0.88"
     * - Integer percentages/bootstrap: "98", "87", "75"
     * 
     * We need to distinguish these from legitimate taxon names. The heuristic is:
     * - If it's a decimal between 0-1, it's a support value
     * - If it's an integer between 50-100, it's likely a bootstrap support value
     * - Otherwise, treat as taxon name (allows "1", "2", "11" as taxa)
     */
    // private boolean isSupportValue(String str, int position, String newickLine) {
    //     try {
    //         double value = Double.parseDouble(str);
            
    //         // Floating point between 0 and 1 is definitely a support value
    //         if(str.contains(".") && (value >= 0 && value <= 1)) {
    //             return true;
    //         }
            
    //         // For integers, use context to determine if it's a support value
    //         // Support values typically appear after ')' or before ':'
    //         if(!str.contains(".") && value >= 50 && value <= 100) {
    //             // Check if this number appears right after ')' (internal node support)
    //             // or right before ':' (support before branch length)
    //             if(position > 0 && newickLine.charAt(position - str.length() - 1) == ')') {
    //                 return true;
    //             }
    //             // Check if followed by ':' (support:length pattern)
    //             if(position < newickLine.length() && newickLine.charAt(position) == ':') {
    //                 return true;
    //             }
    //         }
            
    //         return false;
    //     } catch(NumberFormatException e) {
    //         return false;
    //     }
    // }




    /**
     * Extracts taxon names from a single Newick tree string with support values and branch lengths.
     * - Treat any numeric token that appears immediately AFTER ')' as an internal-node support label.
     * - Otherwise, tokens are taxa (including purely numeric leaf names like "1", "11").
     * - Support labels may be integers, decimals, optional percent (e.g., 98, 0.95, 98%).
     * - Branch lengths after ':' may be in scientific notation.
     */
    private void parseTaxa(String newickLine, Set<String> taxaSet) {
        // Remove whitespace (if you expect quoted labels with spaces, handle quotes instead of stripping)
        newickLine = newickLine.replaceAll("\\s+", "");
        final int n = newickLine.length();
        int i = 0;

        while (i < n) {
            char c = newickLine.charAt(i);

            if (c == '(' || c == ')' || c == ',' || c == ';') {
                i++;
                continue;
            }

            if (c == ':') {
                // Skip branch length (supports scientific notation)
                i++;
                while (i < n) {
                    char bc = newickLine.charAt(i);
                    if (Character.isDigit(bc) || bc == '.' || bc == 'e' || bc == 'E' || bc == '+' || bc == '-') {
                        i++;
                    } else {
                        break;
                    }
                }
                continue;
            }

            // Handle quoted labels: 'A B C' – treat as a single taxon, even if numeric-looking
            if (c == '\'') {
                int start = ++i;
                StringBuilder token = new StringBuilder();
                while (i < n && newickLine.charAt(i) != '\'') {
                    token.append(newickLine.charAt(i));
                    i++;
                }
                // skip closing quote if present
                if (i < n && newickLine.charAt(i) == '\'') i++;
                String label = token.toString();
                if (!label.isEmpty()) {
                    taxaSet.add(label);
                }
                continue;
            }

            // General token (until structural char or ':')
            int start = i;
            StringBuilder token = new StringBuilder();
            while (i < n) {
                char cj = newickLine.charAt(i);
                if (cj == ':' || cj == ')' || cj == ',' || cj == ';' || cj == '(') break;
                token.append(cj);
                i++;
            }
            String label = token.toString();
            if (label.isEmpty()) continue;

            // Context chars around token
            char prev = prevNonSpace(newickLine, start - 1);
            char next = (i < n ? newickLine.charAt(i) : '\0');

            // Decide: support vs taxon
            if (isNumericOrPercent(label) && prev == ')') {
                // numeric token right after ')' = internal-node support
                // do nothing (skip)
            } else {
                // treat as taxon
                taxaSet.add(label);
            }
            // loop continues; i currently at a delimiter or end
        }
    }

    private char prevNonSpace(String s, int idx) {
        while (idx >= 0) {
            char c = s.charAt(idx);
            // we already removed whitespace, but keep this robust
            if (!Character.isWhitespace(c)) return c;
            idx--;
        }
        return '\0';
    }

    /** Accepts integers/decimals with optional % (e.g., 98, 0.95, 98%, 0.95%). */
    private boolean isNumericOrPercent(String s) {
        // Quick path
        if (s.isEmpty()) return false;
        // Optional trailing %
        String core = s.endsWith("%") ? s.substring(0, s.length() - 1) : s;
        if (core.isEmpty()) return false;
        // Match integer/decimal (no exponent for support labels; add if you need it)
        // ^[+-]?(\d+(\.\d+)?|\.\d+)$
        int len = core.length();
        int i = 0;
        if (core.charAt(0) == '+' || core.charAt(0) == '-') {
            if (len == 1) return false;
            i = 1;
        }
        boolean dotSeen = false, digitSeen = false;
        for (; i < len; i++) {
            char c = core.charAt(i);
            if (c == '.') {
                if (dotSeen) return false;
                dotSeen = true;
            } else if (Character.isDigit(c)) {
                digitSeen = true;
            } else {
                return false;
            }
        }
        return digitSeen;
    }








        
    

    /**
     * Reads and standardizes taxon names across all gene trees.
     * 
     * This method performs the first pass through all gene trees to:
     * 1. Collect all unique taxon names that appear in any gene tree
     * 2. Assign consistent IDs to taxa for use throughout the algorithm
     * 3. Create the taxaMap used for consistent tree parsing
     * 
     * The consistent taxon ID assignment is crucial for Algorithm 2 scoring
     * efficiency, as it enables O(1) taxon lookup operations.
     */
    public Map<String, RealTaxon> readTaxaNames() throws FileNotFoundException{
        Set<String> taxaSet = new HashSet<>();

        // First pass: collect all taxon names
        Scanner scanner = new Scanner(new File(this.path));
        while(scanner.hasNextLine()){
            String line = scanner.nextLine();
            if(line.trim().length() == 0) continue;
            parseTaxa(line, taxaSet);
        }
        scanner.close();

        // Create RealTaxon objects with consistent IDs
        this.taxaMap = new HashMap<>();
        for(var x : taxaSet){
            RealTaxon taxon = new RealTaxon(x);
            taxaMap.put(x, taxon);
        }

        return taxaMap;
    }

    /**
     * Reads and preprocesses all gene trees for wQFM-TREE algorithm.
     * 
     * This method performs the main gene tree preprocessing described in Section 2:
     * 1. Parses each gene tree from Newick format using consistent taxon mapping
     * 2. Optionally resolves polytomies to ensure binary tree structure
     * 3. Calculates tripartition frequencies for Algorithm 1 (consensus tree construction)
     * 4. Prepares trees for efficient Algorithm 2 scoring operations
     * 
     * The preprocessing ensures that:
     * - All trees use consistent taxon IDs for O(1) lookup in scoring
     * - Tree structure supports efficient quartet evaluation
     * - Tripartition data enables consensus-based initial bipartition generation
     * 
     * @param distanceMatrix Optional distance matrix for polytomy resolution
     */
    public void readGeneTrees(double[][] distanceMatrix) throws FileNotFoundException{
        
        int internalNodesCount = 0;

        Scanner scanner = new Scanner(new File(path));

        // Process each gene tree
        while (scanner.hasNextLine()) {

            String line = scanner.nextLine();
            // System.out.println(line);
            if(line.trim().length() == 0) continue;
            
            // Parse tree with consistent taxon mapping
            var tree = new Tree(line, this.taxaMap);
            
            // Optional polytomy resolution for binary tree structure
            if(Config.RESOLVE_POLYTOMY){
                tree.resolveNonBinary(distanceMatrix);
            }

            // System.out.println(tree.root);
            
            // if(tree.checkIfNonBinary()){
            //     continue;
            // }

            // Calculate tripartition frequencies for Algorithm 1
            tree.calculateFrequencies(triPartitions);
            geneTrees.add(tree);
            
            // Track internal node count for complexity analysis
            internalNodesCount += tree.nodes.size() - tree.leavesCount;

        }
        
        scanner.close();
        
        // Build efficient lookup arrays for Algorithm 2 operations
        this.taxonIdToLabel = new String[this.taxaMap.size()];
        this.taxa = new RealTaxon[this.taxaMap.size()];
        this.realTaxaCount = this.taxaMap.size();

        // Populate lookup arrays indexed by taxon ID
        for(var x : this.taxaMap.entrySet()){
            taxonIdToLabel[x.getValue().id] = x.getKey();
            taxa[x.getValue().id] = x.getValue();
        }

        // Output preprocessing statistics
        System.out.println( "taxon count : " + this.taxaMap.size());
        System.out.println("Gene trees count : " + geneTrees.size());
        System.out.println( "total internal nodes : " + internalNodesCount);
        System.out.println( "unique partitions : " + triPartitions.size());

        // if(internalNodesCount == 50000){
        //     System.out.println("No polytomy, skipping");
        //     System.exit(-1);
        // }

    }

    /**
     * Constructor: Initialize gene trees preprocessing from file path.
     * 
     * This constructor sets up the GeneTrees object for processing gene trees
     * from the specified file. The actual preprocessing happens in subsequent
     * method calls to readTaxaNames() and readGeneTrees().
     */
    public GeneTrees(String path) throws FileNotFoundException{

        this.geneTrees = new ArrayList<>();
        this.triPartitions = new HashMap<>();
        this.path = path;
    }


    /**
     * Constructor: Initialize with existing taxon mapping.
     * 
     * This constructor is used when taxon mapping is already established,
     * immediately proceeding to gene tree preprocessing.
     */
    public GeneTrees(String path, Map<String, RealTaxon> taxaMap) throws FileNotFoundException{

        this.geneTrees = new ArrayList<>();
        this.triPartitions = new HashMap<>();
        this.path = path;
        this.taxaMap = taxaMap;

        this.readGeneTrees(null);
    }
    
} 


