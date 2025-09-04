import java.util.*;
import src.Tree.Tree;
import src.Tree.TreeNode;
import src.Taxon.RealTaxon;

/**
 * Interactive Quartet Debugger
 * 
 * This tool allows interactive testing of quartet topology detection.
 * - Contains a hardcoded Newick tree (easily editable)
 * - Prompts for 4 taxa names interactively
 * - Shows quartet topology and weight calculation
 * - Runs continuously until user exits
 */
public class QuartetDebugger {
    
    // EDIT THIS NEWICK STRING TO TEST DIFFERENT TREES
    private static final String NEWICK_TREE = "((3,((11,(10,((5,6),(9,(7,8))))),4)),2,1);";
    
    private Tree tree;
    private Map<String, RealTaxon> taxaMap;
    private Scanner scanner;
    
    public QuartetDebugger() {
        this.scanner = new Scanner(System.in);
        initializeTree();
    }
    
    /**
     * Initialize the tree from the hardcoded Newick string
     */
    private void initializeTree() {
        // Extract taxa using the same approach as GeneTrees
        Set<String> taxaSet = new HashSet<>();
        parseTaxa(NEWICK_TREE, taxaSet);
        
        // Create taxa map with consistent IDs
        taxaMap = new HashMap<>();
        for(String taxonName : taxaSet) {
            RealTaxon taxon = new RealTaxon(taxonName); // Use single-arg constructor like GeneTrees
            taxaMap.put(taxonName, taxon);
        }
        
        // Create the tree
        tree = new Tree(NEWICK_TREE, taxaMap);

        tree.setAllBranchLengths(1.0);
        tree.setAllSupportValues(0.9);
        
        System.out.println("=== QUARTET DEBUGGER ===");
        System.out.println("Loaded tree: " + NEWICK_TREE);
        System.out.println("Processed tree: " + tree.getNewickFormat());
        System.out.println("Available taxa: " + taxaSet);
        System.out.println("Total taxa count: " + taxaSet.size());
        System.out.println("Total nodes: " + tree.nodes.size());
        System.out.println("LCA table size: " + tree.nodes.size() + "x" + tree.nodes.size() + " = " + (tree.nodes.size() * tree.nodes.size()) + " entries");
        System.out.println();
    }
    
    /**
     * Extracts taxon names from a single Newick tree string with support values and branch lengths.
     * Uses the same logic as GeneTrees.parseTaxa()
     */
    private void parseTaxa(String newickLine, Set<String> taxaSet) {
        // Remove whitespace
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

            // Handle quoted labels
            if (c == '\'') {
                StringBuilder token = new StringBuilder();
                i++; // Skip opening quote
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

            // Decide: support vs taxon
            if (isNumericOrPercent(label) && prev == ')') {
                // numeric token right after ')' = internal-node support
                // do nothing (skip)
            } else {
                // treat as taxon
                taxaSet.add(label);
            }
        }
    }

    private char prevNonSpace(String s, int idx) {
        while (idx >= 0) {
            char c = s.charAt(idx);
            if (!Character.isWhitespace(c)) return c;
            idx--;
        }
        return '\0';
    }

    /** Accepts integers/decimals with optional % */
    private boolean isNumericOrPercent(String s) {
        if (s.isEmpty()) return false;
        // Optional trailing %
        String core = s.endsWith("%") ? s.substring(0, s.length() - 1) : s;
        if (core.isEmpty()) return false;
        
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
     * Interactive loop for quartet analysis
     */
    public void runInteractiveMode() {
        System.out.println("=== INTERACTIVE QUARTET ANALYSIS ===");
        System.out.println("Enter 4 taxa names separated by spaces (or 'quit' to exit)");
        System.out.println("Example: 1 2 3 4");
        System.out.println();
        
        while(true) {
            System.out.print("Enter 4 taxa: ");
            String input = scanner.nextLine().trim();
            
            if(input.equalsIgnoreCase("quit") || input.equalsIgnoreCase("exit")) {
                System.out.println("Goodbye!");
                break;
            }
            
            String[] taxaNames = input.split("\\s+");
            
            if(taxaNames.length != 4) {
                System.out.println("❌ Error: Please enter exactly 4 taxa names separated by spaces.");
                continue;
            }
            
            // Validate taxa exist in tree
            boolean allValid = true;
            int[] taxaIds = new int[4];
            for(int i = 0; i < 4; i++) {
                RealTaxon taxon = taxaMap.get(taxaNames[i]);
                if(taxon == null) {
                    System.out.println("❌ Error: Taxon '" + taxaNames[i] + "' not found in tree.");
                    System.out.println("   Available taxa: " + taxaMap.keySet());
                    allValid = false;
                    break;
                }
                taxaIds[i] = taxon.id;
            }
            
            if(!allValid) continue;
            
            // Analyze quartet
            analyzeQuartet(taxaNames, taxaIds);
            System.out.println();
        }
    }
    
    /**
     * Analyze a specific quartet and display results
     */
    private void analyzeQuartet(String[] taxaNames, int[] taxaIds) {
        System.out.println("--- Analyzing quartet: " + Arrays.toString(taxaNames) + " ---");
        
        // Try all 3 possible quartet topologies
        String[][] topologies = {
            {taxaNames[0], taxaNames[1], taxaNames[2], taxaNames[3]}, // (ab|cd)
            {taxaNames[0], taxaNames[2], taxaNames[1], taxaNames[3]}, // (ac|bd)
            {taxaNames[0], taxaNames[3], taxaNames[1], taxaNames[2]}  // (ad|bc)
        };
        
        int[][] topologyIds = {
            {taxaIds[0], taxaIds[1], taxaIds[2], taxaIds[3]}, // (ab|cd)
            {taxaIds[0], taxaIds[2], taxaIds[1], taxaIds[3]}, // (ac|bd)
            {taxaIds[0], taxaIds[3], taxaIds[1], taxaIds[2]}  // (ad|bc)
        };
        
        boolean foundValidQuartet = false;
        
        for(int i = 0; i < 3; i++) {
            Tree.QuartetTopology topology = tree.detectQuartetTopology(
                topologyIds[i][0], topologyIds[i][1], topologyIds[i][2], topologyIds[i][3]
            );
            
            if(topology.isValid) {
                foundValidQuartet = true;
                
                // Get the actual pairing from the topology result
                String pair1Taxon1 = getTaxonLabel(topology.pair1[0]);
                String pair1Taxon2 = getTaxonLabel(topology.pair1[1]);
                String pair2Taxon1 = getTaxonLabel(topology.pair2[0]);
                String pair2Taxon2 = getTaxonLabel(topology.pair2[1]);
                
                System.out.println("✅ Valid quartet found!");
                System.out.println("   Topology: ((" + pair1Taxon1 + "," + pair1Taxon2 + "),(" + 
                                 pair2Taxon1 + "," + pair2Taxon2 + "))");
                
                // Calculate and display weight
                double weight = tree.calculateQuartetWeight(topology);
                System.out.println("   Weight: " + String.format("%.6f", weight));
                
                // Show internal nodes info
                System.out.println("   Internal node u (for " + pair1Taxon1 + "," + pair1Taxon2 + 
                                 "): node " + topology.u.index + " at depth " + topology.u.depth + 
                                 ", depthLength " + String.format("%.3f", topology.u.depthLength));
                System.out.println("   Internal node v (for " + pair2Taxon1 + "," + pair2Taxon2 + 
                                 "): node " + topology.v.index + " at depth " + topology.v.depth +
                                 ", depthLength " + String.format("%.3f", topology.v.depthLength));
                
                // Show path lengths (both methods for verification)
                double lengthUA = tree.calculatePathLength(topology.pair1[0], topology.u);
                double lengthUB = tree.calculatePathLength(topology.pair1[1], topology.u);
                double lengthVC = tree.calculatePathLength(topology.pair2[0], topology.v);
                double lengthVD = tree.calculatePathLength(topology.pair2[1], topology.v);
                
                double lengthUA_opt = tree.calculatePathLengthOptimized(topology.pair1[0], topology.u);
                double lengthUB_opt = tree.calculatePathLengthOptimized(topology.pair1[1], topology.u);
                double lengthVC_opt = tree.calculatePathLengthOptimized(topology.pair2[0], topology.v);
                double lengthVD_opt = tree.calculatePathLengthOptimized(topology.pair2[1], topology.v);
                
                System.out.println("   Path lengths (traversal method):");
                System.out.println("     " + pair1Taxon1 + " to u: " + String.format("%.3f", lengthUA));
                System.out.println("     " + pair1Taxon2 + " to u: " + String.format("%.3f", lengthUB));
                System.out.println("     " + pair2Taxon1 + " to v: " + String.format("%.3f", lengthVC));
                System.out.println("     " + pair2Taxon2 + " to v: " + String.format("%.3f", lengthVD));
                
                System.out.println("   Path lengths (optimized formula):");
                System.out.println("     " + pair1Taxon1 + " to u: " + String.format("%.3f", lengthUA_opt));
                System.out.println("     " + pair1Taxon2 + " to u: " + String.format("%.3f", lengthUB_opt));
                System.out.println("     " + pair2Taxon1 + " to v: " + String.format("%.3f", lengthVC_opt));
                System.out.println("     " + pair2Taxon2 + " to v: " + String.format("%.3f", lengthVD_opt));
                
                // Show support product (both methods for verification)
                double supportProduct = tree.calculateSupportProduct(topology.u, topology.v);
                double supportProductOpt = tree.calculateSupportProductOptimized(topology.u, topology.v);
                System.out.println("   Support product along u-v path (traversal): " + String.format("%.6f", supportProduct));
                System.out.println("   Support product along u-v path (optimized): " + String.format("%.6f", supportProductOpt));
                System.out.println("   u depthSupportProductLog: " + String.format("%.6f", topology.u.depthSupportProductLog));
                System.out.println("   v depthSupportProductLog: " + String.format("%.6f", topology.v.depthSupportProductLog));
                
                // Test LCA methods for verification
                TreeNode lcaTraditional = tree.findLCA(topology.u, topology.v);
                TreeNode lcaFast = tree.findLCAFast(topology.u, topology.v);
                System.out.println("   LCA of u,v (traditional): node " + lcaTraditional.index + " at depth " + lcaTraditional.depth);
                System.out.println("   LCA of u,v (O(1) table): node " + lcaFast.index + " at depth " + lcaFast.depth);
                
                break; // Found valid topology, no need to check others
            }
        }
        
        if(!foundValidQuartet) {
            System.out.println("❌ No valid quartet topology found (likely due to polytomy)");
            
            // Show LCA analysis for debugging
            System.out.println("   LCA Analysis:");
            for(int i = 0; i < 4; i++) {
                String refTaxon = taxaNames[i];
                System.out.println("   Reference taxon " + refTaxon + ":");
                for(int j = 0; j < 4; j++) {
                    if(i != j) {
                                                 TreeNode lca = tree.findLCA(taxaIds[i], taxaIds[j]);
                        System.out.println("     LCA with " + taxaNames[j] + ": node " + 
                                         lca.index + " at depth " + lca.depth);
                    }
                }
            }
        }
    }
    
    /**
     * Helper method to get taxon label by ID
     */
    private String getTaxonLabel(int taxonId) {
        for(RealTaxon taxon : taxaMap.values()) {
            if(taxon.id == taxonId) {
                return taxon.label;
            }
        }
        return String.valueOf(taxonId);
    }
    
    /**
     * Main method
     */
    public static void main(String[] args) {
        QuartetDebugger debugger = new QuartetDebugger();
        debugger.runInteractiveMode();
    }
} 