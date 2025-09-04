package src.Tree;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;
import java.util.Stack;
import java.util.List;

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
    
    // LCA table for O(1) LCA queries
    private TreeNode[][] lcaTable;                  // Precomputed LCA table for all node pairs
    private boolean lcaTableBuilt = false;          // Flag to track if LCA table is built


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
        
        // Calculate depths for efficient LCA operations
        calculateDepths();
        
        // Build LCA table for O(1) LCA queries
        buildLCATable();
        
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
     * Calculates depth, cumulative branch lengths, and cumulative support products for all nodes.
     * Depth is measured as number of edges from root (root has depth 0).
     * DepthLength is cumulative branch length from root (root has depthLength 0).
     * DepthSupportProductLog is cumulative log(1-support) from root for efficient support product calculations.
     */
    private void calculateDepths() {
        calculateDepthsUtil(root, 0, 0.0, 0.0);
    }
    
    private void calculateDepthsUtil(TreeNode node, int depth, double depthLength, double depthSupportProductLog) {
        node.setDepth(depth);
        node.setDepthLength(depthLength);
        node.setDepthSupportProductLog(depthSupportProductLog);
        
        if(node.childs != null) {
            for(TreeNode child : node.childs) {
                // Calculate child's cumulative values
                double childDepthLength = depthLength + child.branchLength;
                
                // Calculate child's cumulative log support product
                double childSupportProductLog;
                double oneMinusSupport = 1.0 - child.support;
                if(oneMinusSupport <= 0.0) {
                    // Handle log(0) case - set to very negative value
                    childSupportProductLog = depthSupportProductLog - 9999.0;
                } else {
                    childSupportProductLog = depthSupportProductLog + Math.log(oneMinusSupport);
                }
                
                calculateDepthsUtil(child, depth + 1, childDepthLength, childSupportProductLog);
            }
        }
    }
    
    /**
     * Finds the Lowest Common Ancestor (LCA) of two nodes.
     * Uses the depth information to efficiently find LCA by moving up from deeper node.
     */
    public TreeNode findLCA(TreeNode node1, TreeNode node2) {
        if(node1 == null || node2 == null) return null;
        
        // Make node1 the deeper node
        if(node1.depth < node2.depth) {
            TreeNode temp = node1;
            node1 = node2;
            node2 = temp;
        }
        
        // Move node1 up until both are at same depth
        while(node1.depth > node2.depth) {
            node1 = node1.parent;
        }
        
        // Move both up until they meet
        while(node1 != node2) {
            node1 = node1.parent;
            node2 = node2.parent;
        }
        
        return node1;
    }
    
    /**
     * Finds LCA using taxon IDs (convenience method).
     */
    public TreeNode findLCA(int taxonId1, int taxonId2) {
        if(!isTaxonPresent(taxonId1) || !isTaxonPresent(taxonId2)) {
            return null;
        }
        return findLCA(leaves[taxonId1], leaves[taxonId2]);
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
        
        // Rebuild LCA table after tree structure changes
        buildLCATable();
        calculateDepths();
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

        // System.out.println("originalSupport: " + originalSupport);
        // System.out.println("originalLength: " + originalLength);
        
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
        closest.support = 0;              // Full support for the "broken" side
        closest.branchLength = 0;        // Epsilon length for the "broken" side

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

    /**
     * Quartet topology representation
     */
    public static class QuartetTopology {
        public int[] pair1;  // First pair of taxa (a,b)
        public int[] pair2;  // Second pair of taxa (c,d)
        public TreeNode u;   // Internal node where pair1 meets
        public TreeNode v;   // Internal node where pair2 meets
        public boolean isValid; // Whether this forms a valid quartet
        
        public QuartetTopology(int a, int b, int c, int d, TreeNode u, TreeNode v) {
            this.pair1 = new int[]{a, b};
            this.pair2 = new int[]{c, d};
            this.u = u;
            this.v = v;
            this.isValid = true;
        }
        
        public QuartetTopology() {
            this.isValid = false;
        }
    }
    
    /**
     * Detects quartet topology for four given taxa.
     * 
     * For taxa a,b,c,d, determines which pairing ((a,b),(c,d)) vs ((a,c),(b,d)) vs ((a,d),(b,c))
     * forms a valid quartet in this tree. Also finds the internal nodes u and v where pairs meet.
     * 
     * Algorithm:
     * 1. For each taxon, find LCA with the other three
     * 2. Check if the lowest LCA is unique (indicates valid quartet)
     * 3. Find internal nodes u and v for the valid pairing
     */
    public QuartetTopology detectQuartetTopology(int a, int b, int c, int d) {
        // Check if all taxa are present
        if(!isTaxonPresent(a) || !isTaxonPresent(b) || !isTaxonPresent(c) || !isTaxonPresent(d)) {
            return new QuartetTopology(); // Invalid quartet
        }
        
        int[] taxa = {a, b, c, d};
        
        // Try each taxon as the reference to find unique lowest LCA
        for(int i = 0; i < 4; i++) {
            int ref = taxa[i];
            int[] others = new int[3];
            int idx = 0;
            for(int j = 0; j < 4; j++) {
                if(j != i) others[idx++] = taxa[j];
            }
            
            // Find LCAs of reference with the other three (using O(1) fast LCA)
            TreeNode lca1 = findLCAFast(ref, others[0]);
            TreeNode lca2 = findLCAFast(ref, others[1]);
            TreeNode lca3 = findLCAFast(ref, others[2]);
            
            // Check if one LCA is uniquely deepest (larger depth = further from root, closer to leaves)
            TreeNode[] lcas = {lca1, lca2, lca3};
            int maxDepth = Math.max(Math.max(lca1.depth, lca2.depth), lca3.depth);
            
            int countAtMaxDepth = 0;
            int maxDepthIndex = -1;
            for(int k = 0; k < 3; k++) {
                if(lcas[k].depth == maxDepth) {
                    countAtMaxDepth++;
                    maxDepthIndex = k;
                }
            }
            
            // If exactly one LCA is at maximum depth, we found our quartet
            if(countAtMaxDepth == 1) {
                int partner = others[maxDepthIndex];
                int[] remaining = new int[2];
                idx = 0;
                for(int k = 0; k < 3; k++) {
                    if(k != maxDepthIndex) remaining[idx++] = others[k];
                }
                
                // We have pairs: (ref, partner) and (remaining[0], remaining[1])
                // this v is for ref and partner 
                TreeNode v = lcas[maxDepthIndex]; // LCA of the pair with deepest depth
                
                // Find u carefully
                TreeNode u = findInternalNodeU(remaining[0], remaining[1], ref, partner, v);
                
                return new QuartetTopology(ref, partner, remaining[0], remaining[1], v, u);
            }
        }
        
        return new QuartetTopology(); // No valid quartet found (polytomy)
    }
    
    /**
     * Finds internal node u for quartet topology.
     * Given that v is LCA of one pair, finds u for the other pair.
     */
    private TreeNode findInternalNodeU(int a, int b, int c, int d, TreeNode v) {
        TreeNode lcaAC = findLCAFast(a, c);
        TreeNode lcaBC = findLCAFast(b, c);
        
        if(lcaAC != lcaBC) {
            // Return the one with smaller depth (closer to leaves)
            return (lcaAC.depth > lcaBC.depth) ? lcaAC : lcaBC;
        } else {
            // Both are same, so u is LCA of the first pair
            return findLCAFast(a, b);
        }
    }

    /**
     * Calculates the sum of branch lengths between any two nodes in the tree.
     * Uses LCA traversal approach (kept for compatibility/verification).
     */
    public double calculatePathLength(TreeNode node1, TreeNode node2) {
        if(node1 == null || node2 == null) return 0.0;
        if(node1 == node2) return 0.0; // Same node, no distance
        
        TreeNode lca = findLCA(node1, node2);
        double totalLength = 0.0;
        
        // Path from node1 to LCA
        TreeNode current = node1;
        while(current != null && current != lca) {
            totalLength += current.branchLength;
            current = current.parent;
        }
        
        // Path from node2 to LCA
        current = node2;
        while(current != null && current != lca) {
            totalLength += current.branchLength;
            current = current.parent;
        }
        
        return totalLength;
    }
    
    /**
     * Calculates the sum of branch lengths from a taxon (leaf) to a given node.
     * Convenience method for taxon ID to node calculations.
     */
    public double calculatePathLength(int taxonId, TreeNode targetNode) {
        if(!isTaxonPresent(taxonId)) return 0.0;
        return calculatePathLength(leaves[taxonId], targetNode);
    }
    
    /**
     * Calculates the product of (1 - support) values along the path between two nodes.
     * This is used in the quartet weight formula: product of (1-s(e)) for all branches e along path u,v.
     */
    public double calculateSupportProduct(TreeNode node1, TreeNode node2) {
        if(node1 == null || node2 == null) return 1.0;
        
        // Find path from node1 to node2
        TreeNode lca = findLCA(node1, node2);
        
        double product = 1.0;
        
        // Path from node1 to LCA (excluding LCA)
        TreeNode current = node1;
        while(current != null && current != lca) {
            product *= (1.0 - current.support);
            current = current.parent;
        }
        
        // Path from node2 to LCA (excluding LCA)
        current = node2;
        while(current != null && current != lca) {
            product *= (1.0 - current.support);
            current = current.parent;
        }
        
        return product;
    }

    /**
     * Calculates quartet weight for a given quartet topology.
     * 
     * Formula: e^-(len(u,a) + len(u,b) + len(v,c) + len(v,d)) * (1 - product of (1-s(e)) for path u,v)
     * 
     * Where:
     * - u is internal node where first pair meets
     * - v is internal node where second pair meets  
     * - len(x,y) is sum of branch lengths from x to y
     * - s(e) is support value of branch e
     */
    public double calculateQuartetWeight(QuartetTopology quartet) {
        if(!quartet.isValid) return 0.0;
        
        // Calculate path lengths from internal nodes to leaves (using optimized method)
        double lengthUA = calculatePathLengthOptimized(quartet.pair1[0], quartet.u);
        double lengthUB = calculatePathLengthOptimized(quartet.pair1[1], quartet.u);
        double lengthVC = calculatePathLengthOptimized(quartet.pair2[0], quartet.v);
        double lengthVD = calculatePathLengthOptimized(quartet.pair2[1], quartet.v);
        
        // Total path length component
        double totalLength = lengthUA + lengthUB + lengthVC + lengthVD;
        
        // Support product along path u to v (using optimized method)
        double supportProduct = calculateSupportProductOptimized(quartet.u, quartet.v);
        
        // Final quartet weight formula
        double weight = Math.exp(-totalLength) * (1.0 - supportProduct);
        
        return weight;
    }

    /**
     * Optimized calculation of branch length distance between two nodes.
     * Uses the formula: dist(u,v) = depthLength[u] + depthLength[v] - 2⋅depthLength[lca(u,v)]
     * With O(1) LCA lookup, this is truly O(1).
     */
    public double calculatePathLengthOptimized(TreeNode node1, TreeNode node2) {
        if(node1 == null || node2 == null) return 0.0;
        if(node1 == node2) return 0.0; // Same node, no distance
        
        TreeNode lca = findLCAFast(node1, node2);
        return node1.depthLength + node2.depthLength - 2.0 * lca.depthLength;
    }
    
    /**
     * Optimized path length calculation from taxon to node.
     */
    public double calculatePathLengthOptimized(int taxonId, TreeNode targetNode) {
        if(!isTaxonPresent(taxonId)) return 0.0;
        return calculatePathLengthOptimized(leaves[taxonId], targetNode);
    }

    /**
     * Optimized calculation of support product between two nodes using logarithmic approach.
     * Uses the formula: log(product) = depthSupportProductLog[u] + depthSupportProductLog[v] - 2⋅depthSupportProductLog[lca(u,v)]
     * Then exponentiates the result, with safeguard for very negative values.
     * With O(1) LCA lookup, this is truly O(1).
     */
    public double calculateSupportProductOptimized(TreeNode node1, TreeNode node2) {
        if(node1 == null || node2 == null) return 1.0;
        if(node1 == node2) return 1.0; // Same node, no path
        
        TreeNode lca = findLCAFast(node1, node2);
        double logProduct = node1.depthSupportProductLog + node2.depthSupportProductLog - 2.0 * lca.depthSupportProductLog;
        
        // Handle very negative values (essentially zero product)
        if(logProduct < -100.0) {
            return 0.0;
        }
        
        return Math.exp(logProduct);
    }

    /**
     * Sets all branch lengths in the tree to the specified value.
     * Useful for standardizing trees or testing scenarios.
     * 
     * @param length The branch length value to set for all branches
     */
    public void setAllBranchLengths(double length) {
        for(TreeNode node : nodes) {
            node.setBranchLength(length);
        }
        calculateDepths();
    }
    
    /**
     * Sets all support values in the tree to the specified value.
     * Useful for standardizing trees or testing scenarios.
     * Note: Values will still go through normalization if > 1.0
     * 
     * @param support The support value to set for all branches (should be in [0,1] or [0,100])
     */
    public void setAllSupportValues(double support) {
        for(TreeNode node : nodes) {
            node.setSupport(support);
        }
        // Re-normalize in case the provided support value needs normalization
        normalizeSupportValues();
        calculateDepths();
    }

    /**
     * Builds LCA table for O(1) LCA queries.
     * Uses O(n²) time and space to precompute all pairwise LCAs.
     * 
     * Algorithm:
     * 1. For each node x, find all nodes in subtree of each child
     * 2. Set LCA for pairs (x,v) where v is in subtree of x  
     * 3. Set LCA for cross pairs between different child subtrees
     */
    public void buildLCATable() {
        int n = nodes.size();
        lcaTable = new TreeNode[n][n];
        
        // Initialize diagonal (each node's LCA with itself)
        for(int i = 0; i < n; i++) {
            lcaTable[i][i] = nodes.get(i);
        }
        
        // Build subtree lists for each node
        List<Integer>[] subtree = new List[n];
        for(int i = 0; i < n; i++) {
            subtree[i] = new ArrayList<>();
        }
        buildSubtrees(root, subtree);
        
        // For every node x, assign pairs whose LCA is x
        for(TreeNode x : nodes) {
            int xIndex = x.index;
            
            if(x.childs != null) {
                // Collect subtree lists for each child
                List<List<Integer>> childSubtrees = new ArrayList<>();
                for(TreeNode child : x.childs) {
                    childSubtrees.add(subtree[child.index]);
                    
                    // Pairs (x, v) with v in subtree(child) -> LCA is x
                    for(int v : subtree[child.index]) {
                        if(v != xIndex) { // Don't overwrite diagonal
                            lcaTable[xIndex][v] = x;
                            lcaTable[v][xIndex] = x;
                        }
                    }
                }
                
                // Cross pairs across different child subtrees -> LCA is x
                for(int i = 0; i < childSubtrees.size(); i++) {
                    for(int j = i + 1; j < childSubtrees.size(); j++) {
                        for(int u : childSubtrees.get(i)) {
                            for(int v : childSubtrees.get(j)) {
                                lcaTable[u][v] = x;
                                lcaTable[v][u] = x;
                            }
                        }
                    }
                }
            }
        }
        
        lcaTableBuilt = true;
    }
    
    /**
     * Recursively builds subtree node lists for each node.
     */
    private void buildSubtrees(TreeNode node, List<Integer>[] subtree) {
        subtree[node.index].add(node.index);
        
        if(node.childs != null) {
            for(TreeNode child : node.childs) {
                buildSubtrees(child, subtree);
                subtree[node.index].addAll(subtree[child.index]);
            }
        }
    }
    
    /**
     * O(1) LCA query using precomputed table.
     * Falls back to O(depth) method if table not built.
     */
    public TreeNode findLCAFast(TreeNode node1, TreeNode node2) {
        if(node1 == null || node2 == null) return null;
        
        if(lcaTableBuilt) {
            return lcaTable[node1.index][node2.index];
        } else {
            // Fallback to original method
            return findLCA(node1, node2);
        }
    }
    
    /**
     * O(1) LCA query using taxon IDs.
     */
    public TreeNode findLCAFast(int taxonId1, int taxonId2) {
        if(!isTaxonPresent(taxonId1) || !isTaxonPresent(taxonId2)) {
            return null;
        }
        return findLCAFast(leaves[taxonId1], leaves[taxonId2]);
    }

}
