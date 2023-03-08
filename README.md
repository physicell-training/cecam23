# CECAM 2023 
This repository includes slides, recordings, source code, and resources for the CECAM 2023 workshop (March 2023). 

## Setup Materials 

### PhysiCell set-up guides

#### MacOS Setup
* Slides: [click here](https://github.com/physicell-training/ws2022/raw/main/setup/PhysiCell_ws2022_macOS_setup.pdf)
* Video: [click here](https://youtu.be/Sq9nfKS5U0E)

#### Windows Setup
* Slides: [click here](https://github.com/physicell-training/ws2022/raw/main/setup/PhysiCell_ws2022_Windows_setup.pdf) 
* Video: [click here](https://youtu.be/hIP4JUrViRA)
  
#### Linux Setup
Linux users can most likely skip setup g++ and Python setup, but they should download [PhysiCell 1.10.4 (or later)](https://github.com/MathCancer/PhysiCell/releases/latest) to test that they can compile and run the sample projects. See either the Windows or MacOS guide for details. 

#### Further help
You can get setup help in our [dedicated PhysiCell Slack workspace](https://join.slack.com/t/physicellcomm-sf93727/shared_invite/zt-qj1av6yd-yVeer8VkQaNDjDz7fF00jA). 

### PhysiMESS set-up guide
With PhysiCell setup you should be ready to run PhysiMESS 

* Source code: [click here](https://github.com/physicell-training/cecam23/tree/main/code/PhysiMESS)
* Setup guide: [click here](https://github.com/physicell-training/cecam23/raw/main/code/PhysiMESS/setup/PhysiMESS_setup_guide.pdf)

You can get setup help in our [CECAM Slack workspace](https://join.slack.com/t/cecamecmworks-klu7155/shared_invite/zt-1p9c2arog-9w65U~b5T4R2N~bB2A6uzw).

### Helpful tutorials 

##### Working with PhysiCell Projects 
* Slides: [click here](https://github.com/physicell-training/ws2022/raw/main/sessions/session_01/slides/PhysiCell_ws2022_Session01.pdf)
* Video: [click here](https://youtu.be/fP7-n_RlITU) 

## Key components in code
If you look in the [code](https://github.com/physicell-training/cecam23/tree/main/code) directory, you'll see: 
* [PhysiCell](https://github.com/physicell-training/cecam23/tree/main/code/PhysiCell) - this is a sneak peak of the upcoming 1.11.0 release, which we'llk use throughout the workshop. 
  * You may want to refer to the [release notes](https://github.com/MathCancer/PhysiCell/blob/dev-paul/README.md)
* [PhysiCell Model Builder](https://github.com/physicell-training/cecam23/tree/main/code/PhysiCell-model-builder) - This is an integrated environment for setting up, running, and visualizing PhysiCell simulations. 
* [PhysiMESS](https://github.com/physicell-training/cecam23/tree/main/code/PhysiMESS) - This is a snapshot of the PhysiMESS extension to PhysiCell, which is tailored to ECM modeling. 

We will run most of our tutorials and sessions within these directories, and we'll rely upon this relative pathing. 

## Other suggested tools. 
The setup guides should already suggest ImageMagick. We also recommend [Microsoft VS Code](https://code.visualstudio.com/) as a good, cross-platform code editor. 



## Sessions 

### Monday

#### Plenary (11.30-12.30)
Detangling complex multicellular systems with agent-based modeling - Paul Macklin 
* Slides: [click here](https://github.com/physicell-training/cecam23/raw/main/slides/CECAM_ecm23_macklin_keynote.pdf)
* Recording: [click here](https://epfl.zoom.us/rec/play/n8KFlrt4hWM6kS1lezhrwCNDK44bFMT_rN6w1Wtn4yH9S0inG6_BndE9swmclsPhAXdV4bfimV5MSa11.BFFJiRvErLc8i19y)

#### Tutorial Session 1 (First PhysiCell model) (12.30-13.30)
An hands-on tutorial that introduces the model builder GUI - Lead by Paul Macklin
* Slides: [click here](https://github.com/physicell-training/cecam23/raw/main/slides/CECAM_ecm23_demo1.pdf) 
* Basic Modelling Workflow Recording: [click here](https://epfl.zoom.us/rec/play/Z-83b0kbh4xncpBI0Lajiqf33zNWnEcaI5MAUYhmMVqVVWdIxK10QgBxvt8KUyO_3opYeLszvd9-pDMj.uKQqMjMNaQsWjPJ6)
* Intermediate Modelling Workflow Recording: [click here](https://epfl.zoom.us/rec/play/zoa-wRMNuGA4lQS9TW6cxOQYySV6CVcBLs8Vfgkzg6RFwQiNVCxHmC1wEN69d8pHJWZUtLkWXD5RGXS8.1l69AwjkObR5rzRE)
* Source code: Use `make load PROJ=demo1` to load the full source. 

#### Tutorial Session 2 (Full modeling workflow) (14.30-16.00)
An introduction to signals, behaviours, functions, and a full modelling example - Lead by Paul Macklin 
* Slides: [click here](https://github.com/physicell-training/cecam23/raw/main/slides/CECAM_ecm23_demo2.pdf)
* Full Modelling Workflow Recording: [click here](https://epfl.zoom.us/rec/play/OIRWvIuY50TIOFHeekhTsEiTOb1t1sr7cDpQcMwmJDlRvFV1qi5YBKIxpIgGNFo8eVXfDh2yqdzfp8pA.Zxn_o36aCW8mRnec)
* Source code: Use `make load PROJ=demo2` to load the full source. 

#### 


### Tuesday

#### Tutorial Session 3 (Modelling the ECM as a substrate) (11.30-12.30)
An introduction to substrates and ECM continuum approaches in PhysiCell - Lead by Robyn Shuttleworth
* Slides: [click here](https://github.com/physicell-training/cecam23/raw/main/slides/Physicell_ModellingECM.pdf)
* Recording [click here](https://epfl.zoom.us/rec/play/u-NhBrhNy-JvLb6RH63hplSA6UIhQg3Y3MXYjUMRxTGbWRl3oIxu2i8VIAuAKmDPMD0Rr-bkOVYT-AiN.LVxaQSwDh8SSHDf4)
* Source code: required files [pyMCDS.py](https://github.com/physicell-training/cecam23/blob/main/code/PhysiCell/pyMCDS.py) and [Simulation_Printing.ipynb](https://github.com/physicell-training/cecam23/blob/main/code/PhysiCell/Simulation_Printing.ipynb)

#### Tutorial Session 4 (Introducing PhysiMESS) (13.30-15.00)
An introduction to PhysiMESS and modelling the ECM components as additional agents - Lead by Cicely Macnamara 

This tutorial will follow on from intial setup as per [setup guide](https://github.com/physicell-training/cecam23/raw/main/code/PhysiMESS/setup/PhysiMESS_setup_guide.pdf)
* Slides [click here](https://github.com/physicell-training/cecam23/raw/main/code/PhysiMESS/setup/PhysiMESS_walkthrough.pdf)
* Recording [click here](https://epfl.zoom.us/rec/play/V94jxIylwj4YPKehppGDRgMVs_i8ZoqnwTn3LM2xgQSn-h_Jg7XWxztwrosGt6iMKkyAfxKSlyeBm2ce.KzuGPCLrt0lC_B2R)
* Source Code [click here](https://github.com/physicell-training/cecam23/tree/main/code/PhysiMESS) 

