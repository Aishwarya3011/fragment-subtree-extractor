# Fragment subtree extractor

## Motivation

Many reconstruction and segmentation pipelines (e.g., CNN/UNet based methods) output large, fragmented structures that can be represented as a single file containing many disconnected or weakly connected rooted components.

In practice, this often results in a giant graph containing thousands of fragments, where meaningful analysis, debugging, or manual proofreading is difficult due to:

* lack of global context
* fragmented connectivity
* high interaction overhead in GUI-based tooling
* difficulty scaling proofreading across large datasets

This repository provides a lightweight way to convert a large fragmented rooted graph into smaller, anchor centric subgraphs for systematic proofreading and evaluation.

## Approach

The script converts a large fragmented rooted graph into subgraph exports centered around user-provided anchor points.

* Anchor points are read from a CSV (e.g., a list of key locations or manually annotated targets).
* The input graph is assumed to be a rooted forest, where each root is either:
    * a node with parent = -1, or
    * an orphan node whose parent is missing from the file.
* For each root subtree, if any node in that subtree falls within an axis aligned cube centered at an anchor point, the entire subtree is exported as a standalone file for that anchor.

The cube half-width (cube_half) is an adjustable threshold:

* larger values => more comprehensive coverage per anchor
* smaller values => tighter, more selective exports

This approach reduces repeated manual proofreading work by exporting anchor centric fragments that can be rapidly proofread, visualized, and evaluated in external tools.

## Inputs

* Raw fragmented graph file (SWC format) (e.g. ../giant.swc)
* CSV containing anchor coordinates (e.g. ../anchors.csv). File requires x,y,z as header followed by anchor coordinates.

## Output

The script creates:

* One output directory per anchor of the form:
`{x_coordinate}x-{y_coordinate}y-{z_coordinate}z.swc`
where `(x_coordinate, y_coordinate, z_coordinate)` are the anchor coordinates (rounded to integers).
* Inside each anchor directory: multiple SWCs of the form:
`root_########.swc`
Each `root_########.swc` filename stores the original root node ID from the input giant.swc. Inside each exported SWC, nodes are renumbered starting at 1 for SWC validity and standalone viewing.

These SWCs represent fragment subtrees, not complete neuron reconstructions.

### Example output structure

```text
out_dir/
├── 12034x-5542y-880z.swc/
│   ├── root_00000012.swc
│   ├── root_00000487.swc
│   └── root_00001933.swc
├── 12101x-5601y-901z.swc/
│   ├── root_00000003.swc
│   └── root_00001044.swc
└── 13005x-5900y-950z.swc/
    ├── root_00000200.swc
    ├── root_00000201.swc
    └── root_00009012.swc
```

## Limitations and future work

1. Anchor coordinates are required beforehand. The anchor list does not need to be comprehensive and can be updated at any stage. As new anchors are added, new exports can be generated from the same raw graph.
2. This tool does not attempt to solve global fragment connectivity. Connecting fragments into a complete structure requires domain context and/or stronger structural priors. Heuristic merging strategies can introduce new connectivity errors on top of existing reconstruction artifacts.
3. The current implementation is coordinate anchor centric, but the same approach could be extended to other selection strategies (e.g., nearest-root assignment, kNN-based subgraph extraction, radius-based extraction, or feature-based filtering).

## Usage

### 1. Create a python venv

Create a python venv and install any missing dependencies using pip.

### 2. Run the command from the parent directory

```
python raw2swc.py \
--giant_swc ../path/to/giant.swc \
--anchors_csv ../path/to/anchors.csv \
--out_dir ../path/to/save/at \
--cube_half 300
```

300u is a good default for comprehensive subtree coverage but may be adjusted according to image size, brain area, population density etc.
