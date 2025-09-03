package src;

import java.util.*;
import java.io.*;
import src.Tree.Tree;
import src.Taxon.RealTaxon;
import src.PreProcessing.GeneTrees;

/**
 * QuartetAnalyzer: Handles quartet analysis across multiple gene trees.
 * 
 * This class implements quartet-based phylogenetic analysis by:
 * 1. Enumerating all possible quartets from the taxa set
 * 2. Analyzing quartet topology in each gene tree
 * 3. Calculating quartet weights based on branch lengths and support values
 * 4. Aggregating quartet weights across all gene trees
 * 5. Outputting results in the specified format
 */
public class QuartetAnalyzer {
    
    private Map<String, RealTaxon> taxaMap;
    private List<Tree> geneTrees;
    private Map<String, Double> quartetWeights;
    
    public QuartetAnalyzer(Map<String, RealTaxon> taxaMap) {
        this.taxaMap = taxaMap;
        this.geneTrees = new ArrayList<>();
        this.quartetWeights = new HashMap<>();
    }
    
    /**
     * Main method for command line usage.
     * 
     * Usage: java src.QuartetAnalyzer <gene_trees_file> <output_file>
     */
    public static void main(String[] args) throws IOException {
        if(args.length < 2) {
            System.out.println("Usage: java src.QuartetAnalyzer <gene_trees_file> <output_file>");
            System.out.println("Example: java src.QuartetAnalyzer input/all_gt.tre quartet_weights.txt");
            System.exit(-1);
        }
        
        String inputFilePath = args[0];
        String outputFilePath = args[1];
        
        System.out.println("Analyzing quartets from gene trees: " + inputFilePath);
        System.out.println("Output will be written to: " + outputFilePath);
        System.out.println("==========================================");
        
        // Use existing GeneTrees class to read and parse trees
        GeneTrees trees = new GeneTrees(inputFilePath);
        var taxaMap = trees.readTaxaNames();
        
        System.out.println("Found " + taxaMap.size() + " unique taxa: " + taxaMap.keySet());
        System.out.println();
        
        // Read the gene trees (without distance matrix for polytomy resolution)
        trees.readGeneTrees(null);
        
        System.out.println("=== QUARTET ANALYSIS ===");
        System.out.println("Total gene trees: " + trees.geneTrees.size());
        System.out.println("Total taxa: " + trees.realTaxaCount);
        
        // Initialize quartet analyzer
        QuartetAnalyzer analyzer = new QuartetAnalyzer(taxaMap);
        
        // Add all gene trees
        for(Tree geneTree : trees.geneTrees) {
            analyzer.addGeneTree(geneTree);
        }
        
        // Analyze all quartets
        analyzer.analyzeQuartets();
        
        // Print summary
        analyzer.printSummary();
        
        // Write results to output file
        analyzer.writeQuartetWeights(outputFilePath);
        System.out.println("\nQuartet weights successfully written to: " + outputFilePath);
    }
    
    /**
     * Adds a gene tree for quartet analysis.
     */
    public void addGeneTree(Tree geneTree) {
        geneTrees.add(geneTree);
    }
    
    /**
     * Generates all possible combinations of 4 taxa from the taxa set.
     */
    private List<int[]> generateAllQuartets() {
        List<int[]> quartets = new ArrayList<>();
        int n = taxaMap.size();
        
        // Generate all combinations of 4 taxa
        for(int a = 0; a < n; a++) {
            for(int b = a + 1; b < n; b++) {
                for(int c = b + 1; c < n; c++) {
                    for(int d = c + 1; d < n; d++) {
                        quartets.add(new int[]{a, b, c, d});
                    }
                }
            }
        }
        
        return quartets;
    }
    
    /**
     * Converts quartet topology to a standardized string representation.
     * Format: ((pair1_taxon1,pair1_taxon2),(pair2_taxon1,pair2_taxon2))
     */
    private String quartetToString(Tree.QuartetTopology quartet) {
        if(!quartet.isValid) return null;
        
        // Get taxon labels instead of IDs
        String[] p1Labels = {getTaxonLabel(quartet.pair1[0]), getTaxonLabel(quartet.pair1[1])};
        String[] p2Labels = {getTaxonLabel(quartet.pair2[0]), getTaxonLabel(quartet.pair2[1])};
        
        // Sort pairs lexicographically by label
        Arrays.sort(p1Labels);
        Arrays.sort(p2Labels);
        
        // Ensure lexicographic ordering of pairs
        if(p1Labels[0].compareTo(p2Labels[0]) > 0 || 
           (p1Labels[0].equals(p2Labels[0]) && p1Labels[1].compareTo(p2Labels[1]) > 0)) {
            String[] temp = p1Labels;
            p1Labels = p2Labels;
            p2Labels = temp;
        }
        
        return String.format("((%s,%s),(%s,%s))", p1Labels[0], p1Labels[1], p2Labels[0], p2Labels[1]);
    }
    
    /**
     * Helper method to get taxon label by ID.
     */
    private String getTaxonLabel(int taxonId) {
        for(RealTaxon taxon : taxaMap.values()) {
            if(taxon.id == taxonId) {
                return taxon.label;
            }
        }
        return String.valueOf(taxonId); // Fallback to ID if label not found
    }
    
    /**
     * Analyzes all quartets across all gene trees and accumulates weights.
     */
    public void analyzeQuartets() {
        List<int[]> allQuartets = generateAllQuartets();
        
        System.out.println("Analyzing " + allQuartets.size() + " quartets across " + geneTrees.size() + " gene trees...");
        
        for(int[] quartet : allQuartets) {
            int a = quartet[0], b = quartet[1], c = quartet[2], d = quartet[3];
            
            // Try all 3 possible quartet topologies for this set of 4 taxa
            analyzeQuartetInAllTrees(a, b, c, d);
        }
    }
    
    /**
     * Analyzes a specific quartet (4 taxa) across all gene trees.
     * Tries all 3 possible topologies: (ab|cd), (ac|bd), (ad|bc)
     */
    private void analyzeQuartetInAllTrees(int a, int b, int c, int d) {
        // The 3 possible quartet topologies for taxa a,b,c,d
        int[][] topologies = {
            {a, b, c, d}, // (ab|cd)
            {a, c, b, d}, // (ac|bd) 
            {a, d, b, c}  // (ad|bc)
        };
        
        for(Tree geneTree : geneTrees) {
            // Try to detect which topology exists in this gene tree
            for(int[] topology : topologies) {
                Tree.QuartetTopology detected = geneTree.detectQuartetTopology(
                    topology[0], topology[1], topology[2], topology[3]
                );
                
                if(detected.isValid) {
                    // Calculate weight for this quartet in this gene tree
                    double weight = geneTree.calculateQuartetWeight(detected);
                    
                    // Convert to string representation and accumulate weight
                    String quartetString = quartetToString(detected);
                    if(quartetString != null) {
                        quartetWeights.merge(quartetString, weight, Double::sum);
                    }
                    break; // Found valid topology, no need to try others for this tree
                }
            }
        }
    }
    
    /**
     * Outputs quartet weights in the specified format to a file.
     * Format: ((1,10),(2,3)); 1.234
     */
    public void writeQuartetWeights(String outputFile) throws IOException {
        try (PrintWriter writer = new PrintWriter(new FileWriter(outputFile))) {
            // Sort quartets for consistent output
            List<String> sortedQuartets = new ArrayList<>(quartetWeights.keySet());
            Collections.sort(sortedQuartets);
            
            for(String quartet : sortedQuartets) {
                double weight = quartetWeights.get(quartet);
                if(weight > 0.0) { // Only output quartets with positive weights
                    writer.printf("%s; %.6f%n", quartet, weight);
                }
            }
        }
    }
    
    /**
     * Prints quartet analysis summary.
     */
    public void printSummary() {
        long totalQuartets = quartetWeights.size();
        long nonZeroQuartets = quartetWeights.values().stream()
            .mapToLong(w -> w > 0.0 ? 1 : 0)
            .sum();
        
        System.out.println("Quartet Analysis Summary:");
        System.out.println("Total quartets analyzed: " + totalQuartets);
        System.out.println("Quartets with non-zero weights: " + nonZeroQuartets);
        System.out.println("Average weight: " + 
            quartetWeights.values().stream().mapToDouble(Double::doubleValue).average().orElse(0.0));
    }
} 