import java.io.*;
import java.util.*;

/**
 * NewickCleaner: Efficiently removes branch lengths and support values from Newick format trees.
 * 
 * This utility performs single-pass string processing to clean gene trees without
 * building intermediate tree structures, making it more reliable and faster for
 * large datasets.
 * 
 * Usage: java NewickCleaner <input_file> <output_file>
 * 
 * Handles various Newick formats:
 * - (A:1.0,B:2.0)0.95:3.0; -> (A,B);
 * - ((A:0.1,B:0.2)0.8:0.3,C:0.4)1.0:0.0; -> ((A,B),C);
 * - A:1.5; -> A;
 */
public class NewickCleaner {
    
    /**
     * Cleans a single Newick string by removing all branch lengths and support values.
     * 
     * Algorithm:
     * 1. Use regex-based approach to identify and remove branch info patterns
     * 2. Handle support values after closing parentheses: )0.95:1.0 -> )
     * 3. Handle branch lengths after taxa/nodes: A:1.0 -> A
     * 4. Preserve quoted taxon names and all structural characters
     * 
     * @param newick The input Newick string
     * @return Cleaned Newick string without branch info
     */
    public static String cleanNewickString(String newick) {
        if(newick == null || newick.trim().isEmpty()) {
            return newick;
        }
        
        String cleaned = newick.trim();
        
        // Pattern 1: Remove support values and branch lengths after closing parentheses
        // Matches: )0.95:1.0 or )0.95 or ):1.0 -> )
        // This handles internal node support values and branch lengths
        cleaned = cleaned.replaceAll("\\)[0-9]*\\.?[0-9]*:?[0-9]*\\.?[0-9]*", ")");
        
        // Pattern 2: Remove branch lengths after any content (taxon names, quoted names, etc.)
        // Matches: A:1.0 or 'Species A':1.0 -> A or 'Species A'
        // Use non-greedy matching to avoid consuming too much
        cleaned = cleaned.replaceAll("([^(),;\\s]+|'[^']*'|\"[^\"]*\"):[0-9]*\\.?[0-9]*", "$1");
        
        // Pattern 3: Clean up any remaining standalone support values after closing parentheses
        // This is a more targeted approach to handle leftover support values
        StringBuilder result = new StringBuilder();
        int i = 0;
        int n = cleaned.length();
        
        while(i < n) {
            char curr = cleaned.charAt(i);
            
            if(curr == ')') {
                result.append(curr);
                i++;
                
                // Skip any immediately following support values
                while(i < n && Character.isWhitespace(cleaned.charAt(i))) {
                    i++; // Skip whitespace
                }
                while(i < n && (Character.isDigit(cleaned.charAt(i)) || cleaned.charAt(i) == '.')) {
                    i++; // Skip support value
                }
            }
            else {
                result.append(curr);
                i++;
            }
        }
        
        return result.toString();
    }
    
    /**
     * Processes a file containing multiple Newick trees.
     * Each line is treated as a separate tree.
     * 
     * @param inputFile Path to input file containing Newick trees
     * @param outputFile Path to output file for cleaned trees
     * @throws IOException If file operations fail
     */
    public static void cleanNewickFile(String inputFile, String outputFile) throws IOException {
        try (BufferedReader reader = new BufferedReader(new FileReader(inputFile));
             PrintWriter writer = new PrintWriter(new FileWriter(outputFile))) {
            
            String line;
            int lineNumber = 0;
            
            while((line = reader.readLine()) != null) {
                lineNumber++;
                
                // Skip empty lines and comments
                line = line.trim();
                if(line.isEmpty() || line.startsWith("#")) {
                    writer.println(line);
                    continue;
                }
                
                try {
                    String cleaned = cleanNewickString(line);
                    writer.println(cleaned);
                    
                    // Progress feedback for large files
                    if(lineNumber % 1000 == 0) {
                        System.out.println("Processed " + lineNumber + " trees...");
                    }
                } catch(Exception e) {
                    System.err.println("Error processing line " + lineNumber + ": " + line);
                    System.err.println("Error: " + e.getMessage());
                    // Write original line as fallback
                    writer.println(line);
                }
            }
            
            System.out.println("Successfully processed " + lineNumber + " trees.");
        }
    }
    
    /**
     * Main method - handles command line arguments and file processing.
     * 
     * @param args Command line arguments: [input_file] [output_file]
     */
    public static void main(String[] args) {
        if(args.length != 2) {
            System.err.println("Usage: java NewickCleaner <input_file> <output_file>");
            System.err.println("  input_file:  Path to file containing Newick format gene trees");
            System.err.println("  output_file: Path to output file for cleaned trees");
            System.exit(1);
        }
        
        String inputFile = args[0];
        String outputFile = args[1];
        
        // Validate input file exists
        File input = new File(inputFile);
        if(!input.exists()) {
            System.err.println("Error: Input file does not exist: " + inputFile);
            System.exit(1);
        }
        
        if(!input.canRead()) {
            System.err.println("Error: Cannot read input file: " + inputFile);
            System.exit(1);
        }
        
        // Create output directory if needed
        File output = new File(outputFile);
        File outputDir = output.getParentFile();
        if(outputDir != null && !outputDir.exists()) {
            if(!outputDir.mkdirs()) {
                System.err.println("Error: Cannot create output directory: " + outputDir.getPath());
                System.exit(1);
            }
        }
        
        try {
            System.out.println("Cleaning Newick trees...");
            System.out.println("Input file: " + inputFile);
            System.out.println("Output file: " + outputFile);
            
            long startTime = System.currentTimeMillis();
            cleanNewickFile(inputFile, outputFile);
            long endTime = System.currentTimeMillis();
            
            System.out.println("Cleaning completed in " + (endTime - startTime) + " ms");
            
        } catch(IOException e) {
            System.err.println("Error processing files: " + e.getMessage());
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    /**
     * Test method to validate the cleaner with sample inputs.
     * Can be called independently for testing.
     */
    public static void runTests() {
        System.out.println("Running NewickCleaner tests...");
        
        String[] testCases = {
            "(A:1.0,B:2.0)0.95:3.0;",
            "((A:0.1,B:0.2)0.8:0.3,C:0.4)1.0:0.0;",
            "A:1.5;",
            "(A,B:1.0,C:2.0);",
            "((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6):0.0;",
            "(('Species A':0.1,'Species B':0.2)0.9:0.3,C:0.4);",
            ""
        };
        
        String[] expected = {
            "(A,B);",
            "((A,B),C);",
            "A;",
            "(A,B,C);",
            "((A,B),(C,D));",
            "(('Species A','Species B'),C);",
            ""
        };
        
        boolean allPassed = true;
        for(int i = 0; i < testCases.length; i++) {
            String result = cleanNewickString(testCases[i]);
            if(!result.equals(expected[i])) {
                System.out.println("Test " + (i+1) + " FAILED:");
                System.out.println("  Input:    " + testCases[i]);
                System.out.println("  Expected: " + expected[i]);
                System.out.println("  Got:      " + result);
                allPassed = false;
            } else {
                System.out.println("Test " + (i+1) + " PASSED");
            }
        }
        
        if(allPassed) {
            System.out.println("All tests passed!");
        } else {
            System.out.println("Some tests failed!");
        }
    }
} 