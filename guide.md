# wQFM-TREE Codebase Reading Guide

This guide provides an optimal reading order for understanding the wQFM-TREE codebase, designed to build your knowledge progressively from foundational concepts to complex algorithmic implementations.

## Phase 1: Foundation and Configuration
**Start here to understand the basic structure and setup**

### 1. **Config.java**
- **Why first:** Establishes all algorithmic choices and parameters
- **What to focus on:** 
  - Weight normalization strategies (FLAT vs NESTED)
  - Scoring equation variations 
  - Non-quartet handling types (A vs B)
  - FM algorithm iteration limits
- **Key insight:** This file controls how the paper's algorithms are configured in practice

### 2. **Main.java** 
- **Why next:** Shows the complete algorithm workflow
- **What to focus on:**
  - Four main steps: preprocessing → consensus partition → wQFM-TREE execution → output
  - Command line argument handling
  - High-level algorithm orchestration
- **Key insight:** Provides the "big picture" view of how everything connects

## Phase 2: Data Structures and Preprocessing
**Build understanding of how data flows through the system**

### 3. **src/Taxon/DummyTaxon.java**
- **Why here:** Core concept for divide-and-conquer approach
- **What to focus on:**
  - Tree structure representation of dummy taxa
  - Weight normalization implementation (Section 2.3)
  - Flattened vs hierarchical access patterns
- **Key insight:** Understanding dummy taxa is crucial for Algorithm 2 scoring

### 4. **src/Tree/Tree.java**
- **Why next:** Fundamental data structure for gene trees
- **What to focus on:**
  - Newick parsing and tree construction
  - Polytomy resolution mechanisms
  - Tree balancing and re-rooting operations
- **Key insight:** Gene trees are the primary input for all scoring calculations

### 5. **src/PreProcessing/GeneTrees.java**
- **Why here:** See how raw input becomes usable data
- **What to focus on:**
  - Taxa name extraction and standardization
  - Gene tree parsing with consistent taxon mapping
  - Tripartition frequency calculation for Algorithm 1
- **Key insight:** Data preprocessing enables efficient Algorithm 2 operations

## Phase 3: Partition Management
**Understand how taxa are organized during divide-and-conquer**

### 6. **src/DSPerLevel/TaxaPerLevelWithPartition.java**
- **Why now:** Central state management for subproblems
- **What to focus on:**
  - Taxa partition representation (real vs dummy)
  - Weight coefficient storage and lookup
  - Efficient taxon transfer operations for FM algorithm
- **Key insight:** This class maintains all state needed for Algorithm 2 scoring

## Phase 4: Core Algorithms
**Dive into the main algorithmic implementations**

### 7. **src/InitialPartition/ConsensusTreePartition.java**
- **Why start algorithms here:** Implements Algorithm 1 from paper
- **What to focus on:**
  - Consensus tree construction from gene trees
  - Candidate bipartition evaluation
  - Dummy taxa assignment using weighted sums
- **Key insight:** Eliminates O(n⁴) quartet sorting from original wQFM

### 8. **src/ScoreCalculator/NumSatCalculatorNodeEv2.java**
- **Why next:** Implements Algorithm 2 mathematical formulations
- **What to focus on:**
  - Restructured scoring equation from Section 2.5
  - Inter-branch and single-branch weight calculations
  - Type A vs Type B unresolved quartet strategies
- **Key insight:** Direct scoring from gene trees without quartet enumeration

### 9. **src/DSPerLevel/BookKeepingPerLevelv2.java**
- **Why here:** Coordinates Algorithm 2 across all gene trees
- **What to focus on:**
  - Score aggregation from multiple gene trees
  - Gain calculations for FM algorithm
  - Dummy taxon weight normalization coordination
- **Key insight:** Central orchestrator for efficient scoring operations

## Phase 5: Main Algorithm Flow
**See how everything comes together**

### 10. **src/QFM2.java**
- **Why last:** Main algorithm implementing divide-and-conquer
- **What to focus on:**
  - Recursive bipartition refinement using FM algorithm
  - Divide step with dummy taxa creation
  - Conquer step combining subproblem solutions
  - Iterative improvement with gain calculations
- **Key insight:** Complete implementation of wQFM-TREE's core approach

## Reading Tips for Each Phase

### **Phase 1-2: Building Mental Model**
- Read comments to understand paper connections
- Don't worry about implementation details yet
- Focus on understanding data flow and overall structure

### **Phase 3-4: Algorithm Understanding** 
- Pay attention to mathematical formulations in comments
- Cross-reference with paper sections mentioned in documentation
- Understand how dummy taxa enable efficient computation

### **Phase 5: Integration**
- See how all components work together
- Understand the divide-and-conquer recursion
- Appreciate the efficiency gains over original wQFM

## Key Concepts to Track Throughout

1. **Weight Normalization:** How dummy taxa maintain unit contribution
2. **Quartet Classification:** Satisfied, violated, deferred, unresolved
3. **Algorithm 1 vs Algorithm 2:** Initial partitioning vs direct scoring
4. **FM Algorithm:** Iterative bipartition improvement
5. **Divide-and-Conquer:** How subproblems are created and combined

This reading order ensures you build foundational understanding before tackling complex algorithmic details, making the codebase approachable and setting you up for confident future modifications.
