# PIP5K-Ras-bistability-symmetry-breaking

This repository hosts the code that was developed for the article: _"PIP5K-Ras bistability triggers plasma membrane symmetry breaking to define cellularpolarity and regulate migration"_ (DOI: [10.1101/2024.09.15.613115](https://doi.org/10.1101/2024.09.15.613115)). 

Details on the computational simulations and image processing workflows are available in the article's _Methods_ information.

This repository contains primarily two sections:

- `Computational_Simulations`: Contains the code for the reaction-diffusion based excitable network model and level-set method based viscoelastic model for cell deformation. 

- `Image_Analysis`: Contains the code for image processing workflows, including hierarchical cell classification and cell migration analysis pipeline.

## System Requirements

The code in the `Computational_Simulations` section is implemented in MATLAB 2025a (MathWorks, Natick, MA) inside MacOS. `SDE Toolbox` and `A Toolbox of Level Set Methods` are required for the reaction-diffusion and level-set method based simulations, respectively. 

The code in the `Image_Analysis` section is implemented in Python 3.10 inside Ubuntu 22.04 LTS (except the Jython script for TrackMate, as indicated). Please see the [requirements.txt](Image_Analysis/requirements.txt) file for the list of exact versions of the required Python packages. 

## Installation guide

1. Clone the repository to your local machine using the following command:

   ```bash
   git clone https://github.com/tatsatb/PIP5K-Ras-bistability-symmetry-breaking/
    ```

2. For the `Computational_Simulations` section, ensure that you have MATLAB 2025a installed on your machine. Install the required toolboxes (`SDE Toolbox` and `A Toolbox of Level Set Methods`) from the MathWorks website. Note that the code should be run inside MacOS and may require adjustments for other operating systems.

3. For the `Image_Analysis` section, ensure that you have Python 3.10 installed and a virtual environment is set up. You can create a virtual environment using you favorite environment manager (e.g., `venv`, `conda`, `pyenv`, `mise`, etc.). Then, install the required Python packages using the following command:


   ```bash
   pip install -r Image_Analysis/requirements.txt
   ```

   Notably, these Jupyter Notebooks should be executed using a Jupyter kernel, running Python 3.10 on a Ubuntu 22.04 LTS system as it may require adjustments for other operating systems or Python versions. Also, note that the Jython script for TrackMate should be run inside Fiji (ImageJ). 

As long as MATALB and Python are properly set up with the required dependencies, the setup should be completed within a few minutes. 

## Demo and Instructions for use

Please go through the `README.md` files within the subdirectories of the `Computational_Simulations` directory, which contain detailed instructions on how to run the MATLAB programs.

Please refer to the notes on the Jupyter Notebooks in the `Image_Analysis` section for detailed instructions on how to use the code for image processing and analysis. 

## Authors

The code in this repository was developed by Debojyoti Biswas, Tatsat Banerjee, Pablo A Iglesias, and Parijat Banerjee (Johns Hopkins University, Baltimore, MD, USA).


## Citation/Restrictions

This program is a free and open-source software (please see the [License](#license) information below for details). However, if you are using this code in your work, please cite our work as:

> Yu Deng, Tatsat Banerjee, Satomi Matsuoka, Debojyoti Biswas, Liz A Kurtz, Jane Borleis, Yu Long, Parijat Banerjee, Huiwang Zhan, Dhiman Sankar Pal, Nathan H Roy, Masahiro Ueda, Pablo A Iglesias, Peter N Devreotes. “_PIP5K-Ras bistability triggers plasma membrane symmetry breaking to define cellular polarity and regulate migration_”, bioRxiv, 2025, 2024.09.15.613115 DOI: [10.1101/2024.09.15.613115](https://doi.org/10.1101/2024.09.15.613115). 


## License

Copyright © 2025 Yu Deng, Tatsat Banerjee, Satomi Matsuoka, Debojyoti Biswas, Liz A Kurtz, Jane Borleis, Yu Long, Parijat Banerjee, Huiwang Zhan, Dhiman Sankar Pal, Nathan H Roy, Masahiro Ueda, Pablo A Iglesias, Peter N Devreotes.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.