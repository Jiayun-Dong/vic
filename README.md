# Variable Importance Could

This repository describes datasets and codes to replicate our experiments in the paper (https://arxiv.org/abs/1901.03209).

## Experiment 1 (Section 5.1 in the paper)

### Data
Our first experiment is based on the COMPAS model and its criminal recidivism data, which can be found at this link: https://corels.eecs.harvard.edu/corels/run.html. For your convenience, we have uploaded `compas.csv` in this repository.

### Code
In our paper, we do variable importance could analysis based on two model classes: logistic regression models and decision tree models.

#### VID for Logistic Regression

The file `compas-logistics.R` and its companion file `ellipsoid.R` contains the code we use for the VIC analysis based on logistic regressions. For more details, please refer to Section 3 (in particular, Section 3.2.1) of our paper.

#### VID for Decision Trees

The file `compas-decisiontrees.R` and its companion file `findCandidates.R` contains the code we use for the VIC analysis based on decision trees. For more details, please refer to Section (in particular, Section 3.2.2) of our paper.

## Experiment 2 (Section 5.2 in the paper)

### Data
Our second experiment uses the in-vehicle coupon recommendation data, which is also studied by Wang et al. (2017) at this link: https://www.researchgate.net/publication/318645502_in-vehicle_coupon_recommendation. The data can also be found there and for your convenience, we have uploaded `coupon.csv` in this repository.

### Code
We do our VIC analysis based on logistic regression models. The code can be found in the file `coupon_analysis.R`.

## Experiment 3 (Section 5.3 in the paper)

### Data
We have collected 1572 images of cats and dogs from ImageNet to do this experiment. We have provided the images here: https://www.dropbox.com/s/6gmai6p2siougea/dogcat-data.zip?dl=0.

### Code
We extract VGG16 features (Simonyan and Zisserman, 2004) and do VIC analysis on the features. The code we used to extract the features can be found in `feature extraction.ipynb`, which is written based on the sample code by François Chollet for the book Deep Learning with Python. 

We save the extracted VGG16 features for the images in `features_small.csv` together with the labels of the images in `outcomes_small.csv`. Our code for VIC analysis can be found in `vgg_vid.R`.

To visualize the selected models (in `model.csv`), we use a sample image and modify it to maximize the likelihood of being identified as a dog/cat by the models. The code for the image modification can be found in `Activate Neuron.ipynb`. The code is written based on the method described in https://blog.keras.io/how-convolutional-neural-networks-see-the-world.html.

