package src;

import java.io.FileNotFoundException;
import src.PreProcessing.GeneTrees;

/**
 * TreeAnalyzer: Simple tool to analyze gene trees using existing wQFM-TREE classes.
 * 
 * This is a lightweight version that reuses the existing GeneTrees class
 * to read and analyze gene tree files without the complex algorithm logic.
 */
public class TreeAnalyzer {
    
    public static void main(String[] args) throws FileNotFoundException {
        if(args.length < 1) {
            System.out.println("Usage: java src.TreeAnalyzer <gene_trees_file>");
            System.out.println("Example: java src.TreeAnalyzer input/gtree_11tax_est_5genes_R1.tre");
            System.exit(-1);
        }
        
        String inputFilePath = args[0];
        System.out.println("Analyzing gene trees from: " + inputFilePath);
        System.out.println("==========================================");
        
        // Use existing GeneTrees class to read and parse trees
        GeneTrees trees = new GeneTrees(inputFilePath);
        var taxaMap = trees.readTaxaNames();
        
        System.out.println("Found " + taxaMap.size() + " unique taxa: " + taxaMap.keySet());
        System.out.println();
        
        // Read the gene trees (without distance matrix for polytomy resolution)
        trees.readGeneTrees(null);
        
        // Report statistics using the existing data structures
        System.out.println("=== GENE TREE ANALYSIS ===");
        System.out.println("Total gene trees: " + trees.geneTrees.size());
        System.out.println("Total taxa: " + trees.realTaxaCount);
        
                 // Analyze each tree
         for(int i = 0; i < trees.geneTrees.size(); i++) {
             var tree = trees.geneTrees.get(i);
             System.out.println("\nTree " + (i+1) + ":");
             System.out.println("  Newick format: " + tree.getNewickFormat());
             System.out.println("  Total nodes: " + tree.nodes.size());
             System.out.println("  Leaf nodes: " + tree.leavesCount);
             System.out.println("  Internal nodes: " + (tree.nodes.size() - tree.leavesCount));
             System.out.println("  Is binary: " + !tree.checkIfNonBinary());
             
             // Report branch information
             System.out.println("  Branch information:");
             for(var node : tree.nodes) {
                 if(node.parent != null) { // Not root
                     String nodeDesc = node.isLeaf() ? 
                         "Leaf " + node.taxon.label : 
                         "Internal node " + node.index;
                     System.out.println("    " + nodeDesc + 
                         " -> support: " + String.format("%.3f", node.support) +
                         ", length: " + String.format("%.3f", node.branchLength));
                 }
             }
         }
        
        // Overall statistics
        System.out.println("\n=== OVERALL STATISTICS ===");
        System.out.println("Total gene trees: " + trees.geneTrees.size());
        System.out.println("Total number of taxa across all trees: " + taxaMap.size());
        System.out.println("Taxa across all trees: " + taxaMap.keySet());
        
        // Count binary vs non-binary trees
        int binaryCount = 0;
        for(var tree : trees.geneTrees) {
            if(!tree.checkIfNonBinary()) {
                binaryCount++;
            }
        }
        System.out.println("Binary trees: " + binaryCount + "/" + trees.geneTrees.size());
        System.out.println("Non-binary trees: " + (trees.geneTrees.size() - binaryCount) + "/" + trees.geneTrees.size());
    }
} 