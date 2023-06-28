# SpacerPlacer

![SpacerPlacer Logo](doc/figures/logo.pdf)


SpacerPlacer is a powerful software for reconstructing ancestral CRISPR spacer arrays along a given phylogenetic tree.
SpacerPlacer is respecting the chronological timeline of spacer insertions as defined by the order of spacers within the CRISPR arrays and uses this additional information to improve the ancestral reconstruction.


## Features

- **Ancestral Spacer Array Reconstruction:** SpacerPlacer allows you to reconstruct ancestral CRISPR spacer arrays based on a given phylogenetic tree.
- **Identification of Ancestral Spacer Acquisitions, Deletions, Duplications, and Rearrangements:** The tool identifies ancestral duplications and rearrangements within the CRISPR array.
- **Spacer Deletion Rate Estimation:** The software provides an estimation of spacer deletion rates within the CRISPR arrays.
- **Analysis of Jointly Lost Spacers:** SpacerPlacer calculates the average length of blocks of jointly lost spacers in the array.
- **optional: Phylogenetic Tree Reconstruction:** In addition to spacer array reconstruction along a given tree, SpacerPlacer can also reconstruct phylogenetic trees based on the CRISPR spacer data.



## Conference Posters

We presented posters at the conference CRISPR2023, showcasing the features and applications of SpacerPlacer and the related CRISPR-evo-inator tool to predict the direction of transcription of CRISPR arrays. You can find the PDF versions of the posters below:

- [SpacerPlacer Poster](doc/posters/crispr2023_spacerplacer.pdf)
- [CRISPR-evo-inator Poster](doc/posters/crispr2023_CRISPRevoinator.pdf)

Please feel free to explore the posters to learn more about SpacerPlacer and its applications.


## Code Availability

We are currently in the process of preparing the code for SpacerPlacer and will soon be adding it to this repository. We appreciate your patience!

To stay up to date with the latest developments, we recommend bookmarking this repository and checking back periodically for updates.

Alternatively, you can reach out to us via email at axel.fehrenbach@uni-tuebingen.de to express your interest in SpacerPlacer. We will be happy to notify you once the code becomes available.

Thank you for your interest in SpacerPlacer!


## Installation


To install SpacerPlacer, follow these steps:

1. Clone the GitHub repository:

   ```bash
   git clone https://github.com/fbaumdicker/SpacerPlacer.git
   ```
   
2. Navigate to the cloned repository and install the required dependencies:

   ```bash
    cd SpacerPlacer
    pip install -r requirements.txt
   ```
