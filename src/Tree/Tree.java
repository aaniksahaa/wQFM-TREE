package src.Tree;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;
import java.util.Stack;

import src.Taxon.RealTaxon;

/**
 * Tree: Core tree data structure for wQFM-TREE algorithm.
 * 
 * This class implements the phylogenetic tree representation used throughout
 * the wQFM-TREE algorithm. It handles various tree operations essential for
 * the divide-and-conquer approach:
 * 
 * 1. Gene tree parsing from Newick format
 * 2. Consensus tree construction and manipulation
 * 3. Tree balancing and re-rooting operations
 * 4. Support for Algorithm 2 scoring through efficient tree traversal
 * 5. Final species tree output generation
 * 
 * Trees in wQFM-TREE contain real taxa as leaves and support dummy taxa
 * through special node markings. The tree structure enables efficient
 * quartet evaluation and the mathematical formulations from Section 2.5.
 */
public class Tree {
    
    // Core tree structure
    public ArrayList<TreeNode> nodes;               // All nodes in tree (internal + leaves)
    public ArrayList<TreeNode> topSortedNodes;      // Nodes in topological order for traversal
    
    public TreeNode root;                           // Root node of the tree
    public Map<String, RealTaxon> taxaMap;          // Mapping from taxon names to RealTaxon objects
    
    // Leaf access optimization for Algorithm 2 scoring
    public TreeNode[] leaves;                       // Fast access to leaf nodes by taxon ID
    public int leavesCount;                         // Number of leaves in tree
    

    /**
     * Creates a new internal or leaf tree node.
     * 
     * This is the fundamental building block for tree construction, used both
     * during Newick parsing and during the conquer phase when combining
     * subproblem solutions.
     */
    public TreeNode addNode(ArrayList<TreeNode> children, TreeNode parent){

        TreeNode nd = new TreeNode().setIndex(nodes.size()).setChilds(children).setParent(parent);
        nodes.add(nd);
        return nd;
    }


    /**
     * Creates a new internal node with specified children.
     * 
     * Used extensively during tree construction and the conquer phase where
     * subproblem solutions are combined through tree grafting operations.
     */
    public TreeNode addInternalNode(ArrayList<TreeNode> children){
        var nd = addNode(children, null);
        for (var x : children)
            x.setParent(nd);
        return nd;
    }

    /**
     * Creates a new leaf node for a real taxon.
     * 
     * Leaf nodes represent the actual species/taxa being analyzed. During
     * Algorithm 2 scoring, these nodes are the evaluation points for the
     * mathematical formulations from Section 2.5.
     */
    public TreeNode addLeaf(RealTaxon taxon){
        var nd = addNode(null, null).setTaxon(taxon);
        return nd;
    }


    /**
     * Parses a phylogenetic tree from Newick format with support values and branch lengths.
     * 
     * Enhanced parser that handles extended Newick format:
     * - Branch lengths: specified after colon (:1.5)
     * - Support values: specified before colon (0.95:1.5)
     * - For leaves without support: defaults to 1.0
     * 
     * Format examples:
     * - (A:1.0,B:2.0)0.95:3.0  -> internal node with support 0.95 and length 3.0
     * - A:1.5                  -> leaf A with length 1.5 and support 1.0 (default)
     */
    private void parseFromNewick(String newickLine){

        int leavesCount = 0;
        nodes = new ArrayList<>();
        Stack<TreeNode> nodeStack = new Stack<>();
        newickLine = newickLine.replaceAll("\\s", "");
        int n = newickLine.length();
        int i = 0;
        

    
        while(i < n){
            char curr = newickLine.charAt(i);
            if(curr == '('){
                // Start of internal node - push sentinel
                nodeStack.push(null);
                i++;
            }
            else if(curr == ')'){
                // End of internal node - collect children and create internal node
                ArrayList<TreeNode> children = new ArrayList<>();
                while(!nodeStack.isEmpty() && nodeStack.peek() != null){
                    children.add(nodeStack.pop());
                }
                if(!nodeStack.isEmpty())
                    nodeStack.pop(); // Remove sentinel
                
                TreeNode internalNode = addInternalNode(children);
                // Parse support value and branch length for this internal node
                i++; // Move past ')'
                i = parseBranchInfo(newickLine, i, n, internalNode);
                
                nodeStack.push(internalNode);
            }
            else if(curr == ',' || curr == ';'){
                // Separators - skip
                i++;
            }
            else{
                // Parse taxon name and create leaf node
                StringBuilder taxonName = new StringBuilder();
                int j = i;
                
                // Extract taxon name (until we hit :, ), or ,)
                while(j < n){
                    char curr_j = newickLine.charAt(j);
                    if(curr_j == ':' || curr_j == ')' || curr_j == ',' || curr_j == ';'){
                        break;
                    }
                    taxonName.append(curr_j);
                    j++;
                }
                
                // Create leaf node
                RealTaxon taxon = this.taxaMap.get(taxonName.toString());
                TreeNode leafNode = addLeaf(taxon);
                leavesCount++;
                
                // Parse branch length for leaf (leaves default to support = 1.0)
                i = j;
                i = parseBranchInfo(newickLine, i, n, leafNode);
                
                nodeStack.push(leafNode);
            }
        }

        this.leavesCount = leavesCount;
    
        root = nodeStack.lastElement();

        // Ensure binary tree structure for efficient Algorithm 2 operations
        if(root.childs.size() > 2)
            balanceRoot();
        
        // Normalize support values to ensure they are in [0,1] range
        normalizeSupportValues();
        
        // Set up data structures for efficient tree operations
        filterLeaves();
        topSort();

        // System.out.println(this.leavesCount);
        // for( i = 0; i < this.leaves.length; ++i){
        //     if(this.leaves[i] == null){
        //         System.out.println("Error: Taxon " + i + " is not present in tree");
        //         System.exit(-1);
        //     }
        //     if(this.leaves[i].taxon.id != i){
        //         System.out.println("Error: Taxon " + i + " not matching");
        //         System.exit(-1);
        //     }
        // }
        // bringLeafsToFront();

    }

    /**
     * Parses branch information (support value and branch length) from Newick format.
     * 
     * Handles the format: [support]:[length] where support is optional for leaves.
     * Examples:
     * - :1.5 -> support=1.0 (default), length=1.5
     * - 0.95:2.0 -> support=0.95, length=2.0
     * - (no info) -> support=1.0, length=0.0 (defaults)
     * 
     * @param newickLine The complete Newick string
     * @param startPos Current position in the string
     * @param endPos End of string
     * @param node The node to set branch info for
     * @return New position after parsing branch info
     */
    private int parseBranchInfo(String newickLine, int startPos, int endPos, TreeNode node) {
        int i = startPos;
        double support = 1.0;  // Default support
        double branchLength = 0.0;  // Default length
        
        // Check if there's branch information
        if(i < endPos && (newickLine.charAt(i) == ':' || Character.isDigit(newickLine.charAt(i)) || newickLine.charAt(i) == '.')) {
            
            // Case 1: Support value followed by colon and length (support:length)
            if(i < endPos && Character.isDigit(newickLine.charAt(i))) {
                StringBuilder supportStr = new StringBuilder();
                while(i < endPos && newickLine.charAt(i) != ':' && newickLine.charAt(i) != ')' && 
                      newickLine.charAt(i) != ',' && newickLine.charAt(i) != ';') {
                    supportStr.append(newickLine.charAt(i));
                    i++;
                }
                
                if(supportStr.length() > 0) {
                    try {
                        support = Double.parseDouble(supportStr.toString());
                    } catch(NumberFormatException e) {
                        support = 1.0; // Default if parsing fails
                    }
                }
            }
            
            // Case 2: Parse branch length after colon
            if(i < endPos && newickLine.charAt(i) == ':') {
                i++; // Skip the colon
                StringBuilder lengthStr = new StringBuilder();
                while(i < endPos && newickLine.charAt(i) != ')' && newickLine.charAt(i) != ',' && 
                      newickLine.charAt(i) != ';') {
                    lengthStr.append(newickLine.charAt(i));
                    i++;
                }
                
                if(lengthStr.length() > 0) {
                    try {
                        branchLength = Double.parseDouble(lengthStr.toString());
                    } catch(NumberFormatException e) {
                        branchLength = 0.0; // Default if parsing fails
                    }
                }
            }
        }
        
        // Set the parsed values
        node.setSupport(support);
        node.setBranchLength(branchLength);
        
        return i;
    }

    /**
     * Normalizes support values across the tree to ensure they are in the 0-1 range.
     * 
     * This method performs a single pass through all nodes to check support values:
     * 1. If all values are in [0,1], no normalization needed
     * 2. If all values are in [0,100], divide by 100 to normalize to [0,1]
     * 3. If any value exceeds 100, throw an error
     * 
     * This ensures consistent support value interpretation across different input formats.
     */
    private void normalizeSupportValues() {
        boolean hasValuesAboveOne = false;
        boolean hasValuesAbove100 = false;
        double maxSupport = 0.0;
        
        // First pass: analyze the range of support values
        for(TreeNode node : nodes) {
            if(node.support > maxSupport) {
                maxSupport = node.support;
            }
            if(node.support > 1.0) {
                hasValuesAboveOne = true;
            }
            if(node.support > 100.0) {
                hasValuesAbove100 = true;
            }
        }
        
        // Check for invalid values (> 100)
        if(hasValuesAbove100) {
            throw new IllegalArgumentException(
                String.format("Invalid support value found: %.3f. Support values must be in range [0,1] or [0,100].", maxSupport)
            );
        }
        
        // Normalize if values are in 0-100 range
        if(hasValuesAboveOne) {
            System.out.println("Normalizing support values from [0,100] to [0,1] range...");
            for(TreeNode node : nodes) {
                node.support = node.support / 100.0;
            }
        }
    }

    /**
     * Creates fast-access array for leaf nodes indexed by taxon ID.
     * 
     * This optimization is crucial for Algorithm 2 scoring, which needs
     * frequent access to leaf nodes by taxon ID during quartet evaluation.
     * The leaves array enables O(1) lookup instead of O(n) tree traversal.
     */
    private void filterLeaves(){
        this.leaves = new TreeNode[this.taxaMap.size()];
        for(var x : nodes){
            if(x.isLeaf()){
                this.leaves[x.taxon.id] = x;
            }
        }
    }

    /**
     * Recursively resolves non-binary internal nodes using distance matrix.
     * 
     * This method converts non-binary (polytomy) nodes into binary nodes using
     * a distance-based heuristic. This is important for Algorithm 2 scoring
     * which assumes binary tree structure for efficient quartet evaluation.
     * 
     * The resolution strategy selects the child subtree with maximum total
     * distance to all other taxa and separates it from the rest.
     */
    public ArrayList<Integer> resolveNonBinaryUtil(TreeNode node, double[][] distanceMatrix){
        if(node.isLeaf()){
            return new ArrayList<>(Arrays.asList(node.taxon.id));
        }
        
        // Collect taxa reachable from each child subtree
        var reachableFromChilds = new ArrayList<ArrayList<Integer>>();
        for(var x : node.childs){
            reachableFromChilds.add(resolveNonBinaryUtil(x, distanceMatrix));
        }
        var allReachableFromChilds = new ArrayList<Integer>();
        for(var x : reachableFromChilds){
            allReachableFromChilds.addAll(x);
        }

        // Resolve polytomy if more than 2 children
        if(node.childs.size() > 2){
            double mxDist = Double.MIN_VALUE;
            int mxIndex = -1;
            
            // Find child subtree with maximum total distance to all others
            for(int i = 0; i < node.childs.size(); ++i){
                double currDist = 0;
                for(var a : allReachableFromChilds){
                    for(var b : reachableFromChilds.get(i)){
                        currDist += distanceMatrix[a][b];
                    }
                }
                if(currDist > mxDist){
                    mxDist = currDist;
                    mxIndex = i;
                }
            }

            // Separate the maximum distance child from others
            var branchI = node.childs.get(mxIndex);
            node.childs.remove(mxIndex);
            
            // Create new internal node for remaining children
            var newNode = addInternalNode(node.childs);
            newNode.setParent(node);
            node.childs = new ArrayList<>();
            node.childs.add(branchI);
            node.childs.add(newNode);
            
            // Recursively resolve the new internal node
            resolveNonBinaryUtil(newNode, distanceMatrix);

        }
        return allReachableFromChilds;
    }

    public void resolveNonBinary(double[][] distanceMatrix){
        // System.out.println(root.childs.size());
        resolveNonBinaryUtil(root, distanceMatrix);
        topSort();
    }

    
    
    // public Tree(String newickLine){
    //     taxaMap = null;
    //     parseFromNewick(newickLine);
    // }

    public Tree(String newickLine, Map<String, RealTaxon> taxaMap){
        this.taxaMap = taxaMap;
        parseFromNewick(newickLine);
    }

    public Tree(){
        taxaMap = null;
        nodes = new ArrayList<>();
    }


    public int dfs(TreeNode node, ArrayList<Integer> subTreeNodeCount){
        if(node.childs == null){
            subTreeNodeCount.set(node.index, 1);
            return 1;
        }
        int res = 0;
        for(var x : node.childs){
            res += dfs(x, subTreeNodeCount);
        }
        subTreeNodeCount.set(node.index, res + 1);
        return res + 1;
    }

    // private void bringLeafsToFront(){
    //     for(int i = 0; i < nodes.size(); ++i){
    //         if(nodes.get(i).isLeaf()){

    //             var curr = nodes.get(i);
    //             var tmp = nodes.get(curr.taxon.id);
                
    //             nodes.set(curr.taxon.id, curr);
    //             nodes.set(i, tmp);
    //             curr.index = curr.taxon.id;
    //             tmp.index = i;
    //         }
    //     }
    // }
    
    // maybe there is problem here, be CAREFUL 
    // TODO: fix this 
    public void reRootTree(TreeNode newRootNode){
        TreeNode newRootP = newRootNode.parent;
        if(newRootP == null) return;
        
        // Store original branch properties of the new root edge
        double originalSupport = newRootNode.support;
        double originalLength = newRootNode.branchLength;
        
        newRootP.childs.remove(newRootNode);

        TreeNode curr = newRootP;
        TreeNode currP, temp;
        currP = curr.parent;
        
        // During rerooting, we need to reverse parent-child relationships
        // and properly handle branch properties
        while(curr != null && currP != null){
            // Store branch properties before changing relationships
            double tempSupport = currP.support;
            double tempLength = currP.branchLength;
            
            curr.childs.add(currP);
            currP.childs.remove(curr);
            temp = currP;
            currP = currP.parent;
            temp.parent = curr;
            
            // Update branch properties: the branch that was from currP to curr
            // is now from curr to temp, so temp inherits the old properties
            temp.support = tempSupport;
            temp.branchLength = tempLength;
            
            curr = temp;
        }
        
        if(newRootNode.isLeaf())
            newRootNode.childs = new ArrayList<>();
        newRootNode.childs.add(newRootP);
        
        // Handle the new root's branch properties
        // The new root has no parent, so we split the original edge
        newRootNode.support = 1.0;      // Root has full support
        newRootNode.branchLength = 1e-6; // Epsilon length for root
        
        // The former parent gets the remaining branch properties
        newRootP.support = originalSupport;
        newRootP.branchLength = originalLength;
        newRootP.parent = newRootNode;
        
        this.root = newRootNode;
    }

    private void balanceRoot(){

        // System.out.println("newick: " + newickFormatUitl(root));

        int n = nodes.size();
        ArrayList<Integer> subTreeNodeCount = new ArrayList<>(n);
        for(int i = 0; i < n; ++i)
            subTreeNodeCount.add(0);
        dfs(root, subTreeNodeCount);
        
        TreeNode closest = root;
        int diff = n;
        int v;
        for(int i = 0; i < n; ++i){
            v = Math.abs(n/2 - subTreeNodeCount.get(i)); 
            if(v < diff){
                diff = v;
                closest = nodes.get(i);
            }
        }
        
        TreeNode closestP = closest.parent;
        if(closestP == null) return; // Already at root
        
        // Store original branch properties of the edge we're breaking
        double originalSupport = closest.support;
        double originalLength = closest.branchLength;

        System.out.println("originalSupport: " + originalSupport);
        System.out.println("originalLength: " + originalLength);
        
        closestP.childs.remove(closest);

        TreeNode curr = closestP;
        TreeNode currP, temp;
        currP = curr.parent;

        // Store branch properties before changing relationships
        double tempSupport = curr.support;
        double tempLength = curr.branchLength;

        double nextTempSupport, nextTempLength;
        
        // Handle branch properties during the rerooting process
        while(curr != null && currP != null){
            nextTempSupport = currP.support;
            nextTempLength = currP.branchLength;
            
            curr.childs.add(currP);
            currP.childs.remove(curr);
            temp = currP;
            currP = currP.parent;
            temp.parent = curr;
            
            // Update branch properties: the branch that was from currP to curr
            // is now from curr to temp, so temp inherits the old properties
            temp.support = tempSupport;
            temp.branchLength = tempLength;
            
            curr = temp;

            tempSupport = nextTempSupport;
            tempLength = nextTempLength;
        }

        ArrayList<TreeNode> arr = new ArrayList<>();
        arr.add(closest);
        arr.add(closestP);

        root = addInternalNode(arr);
        
        // Handle branch properties for the new root configuration
        // The new root is an artificial node breaking the original edge
        root.support = 1.0;      // Root has full support
        root.branchLength = 0.0; // Root has no parent branch
        
        // The two children of the new root should keep their existing branch properties
        // These were correctly parsed from the Newick string - don't overwrite them!
        // here we set originalSupport to the support of the closest node
        // because if we set 1, then it would make all paths along it to be of weight 1
        closest.support = originalSupport;              // Full support for the "broken" side
        closest.branchLength = 1e-6;        // Epsilon length for the "broken" side

        closestP.support = originalSupport;              // Full support for the "broken" side
        closestP.branchLength = originalLength;        // Epsilon length for the "broken" side
        
        // closestP should keep its existing support and length values that were parsed correctly
        // Don't overwrite them - they contain the correct values from the original Newick string

    
    }
    
    
    
    private String newickFormatUitl(TreeNode node){
        StringBuilder sb = new StringBuilder();
        
        if(node.isLeaf()){
            sb.append(node.taxon.label);
        } else {
            sb.append("(");
            for(int i = 0; i < node.childs.size(); ++i){
                sb.append(newickFormatUitl(node.childs.get(i)));
                if(i != node.childs.size() - 1)
                    sb.append(",");
            }
            sb.append(")");
        }
        
        // Add branch information (support:length) if not root
        if(node.parent != null) {
            // For internal nodes, include support value if it's not 1.0
            if(!node.isLeaf()) {
                sb.append(String.format("%.3f", node.support));
            }

            sb.append(":").append(String.format("%.3f", node.branchLength));

            // // Add branch length if it's not 0.0
            // if(node.branchLength != 0.0) {
            //     sb.append(":").append(String.format("%.3f", node.branchLength));
            // }
        }
        
        return sb.toString();
    }


    public String getNewickFormat(){
        return newickFormatUitl(root) + ";";
    }

    public boolean isTaxonPresent(int id){
        return this.leaves[id] != null;
    }
    

    private ArrayList<Integer> getChildrens(TreeNode node, Map<String, TreeNode> triPartitionsMap){
        if(node.isLeaf()){
            ArrayList<Integer> arr = new ArrayList<>();
            arr.add(node.taxon.id);
            return arr;
        }
        ArrayList<ArrayList<Integer>> reachableFromChilds = new ArrayList<>();
        for(var child : node.childs){
            reachableFromChilds.add(getChildrens(child, triPartitionsMap));
        }

        // flag all elems in left and right partition to find elems of third partition
        boolean[] mark = new boolean[this.taxaMap.size()];

        for(var childTaxa : reachableFromChilds){
            for(var x : childTaxa){
                mark[x] = true;
            }
        }
        
        ArrayList<Integer> arr = new ArrayList<>();
        for(int i = 0; i < this.taxaMap.size(); ++i){
            if(!mark[i] && isTaxonPresent(i)){
                arr.add(i);
            }
            // else if(!isTaxonPresent(i)){
            //     System.out.println("No. " + i + " is not present in tree");
            // }
        }

        reachableFromChilds.add(arr);

        String[] partitionStrings = new String[reachableFromChilds.size()];

        for(var x : reachableFromChilds){
            x.sort((a, b) -> a - b);
        }


        for(int i = 0; i < reachableFromChilds.size(); ++i){
            var sb = new StringBuilder();
            reachableFromChilds.get(i).forEach(s -> sb.append(s + '-') );
            partitionStrings[i] = sb.toString();
        }
        
        Arrays.sort(partitionStrings);
        // String key;
        StringBuilder sb = new StringBuilder();
        for(var x : partitionStrings){
            sb.append(x + '|');
        }

        String key = sb.toString();

        if(triPartitionsMap.containsKey(key)){
            triPartitionsMap.get(key).frequency++;
            node.frequency = 0;
        }
        else{
            node.frequency = 1;
            triPartitionsMap.put(key, node);
        }
        reachableFromChilds.remove(reachableFromChilds.size() - 1);
        
        arr = new ArrayList<>();
        for(var childTaxa : reachableFromChilds){
            arr.addAll(childTaxa);
        }
        return arr;
    }

    public void calculateFrequencies(Map<String, TreeNode> triPartitions){
        getChildrens(root, triPartitions);
    }

    private void topSortUtil(TreeNode node, ArrayList<TreeNode> topSort){
        if(node.isLeaf()){
        }
        else{
            for(var x : node.childs){
                topSortUtil(x, topSort);
            }
        }
        topSort.add(node);

    }

    public void topSort(){
        ArrayList<TreeNode> topSort = new ArrayList<>();
        topSortUtil(root, topSort);
        this.topSortedNodes = topSort;
    }

    public boolean checkIfNonBinary(){
        for(var x : nodes){
            if( x.childs != null && x.childs.size() > 2)
                return true;
        }
        return false;
    }

}
