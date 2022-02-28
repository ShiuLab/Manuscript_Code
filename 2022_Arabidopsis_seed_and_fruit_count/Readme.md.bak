# We use two methods to count Arabidopsis seeds: Faster R-CNN and ImageJ. The Faster R-CNN is also used to count the Arabidopsis fruits (siliques).

  * To use our model, please cite our **New Phytologist** paper [High throughput measurement of plant fitness traits with an object detection method using Faster R-CNN](https://www.biorxiv.org/content/10.1101/2021.07.01.450758v2)

  * All the related scripts are [here](https://github.com/ShiuLab/Manuscript_Code/tree/master/2022_Arabidopsis_seed_and_fruit_count/Scripts_for_Faster_R-CNN). 

  * We don't have instructions for fruit counting models, as it is the same as how we do for the seed counting except that the fruit images are not split. If you are interested in using our final fruit counting model to detect and count Arabidopsis fruits, please check [here](https://github.com/ShiuLab/Manuscript_Code/tree/master/2022_Arabidopsis_seed_and_fruit_count/Fruit_counting_model) for the model, and go to the "A. Apply our seed counting model in windows" or "B. Apply our seed counting model in Linux" for the instruction.

  * If you are interested in build your own models, please go to "C. Build your own seed counting models in Linux" for the instruction.

  * If you want to use your imageJ pipeline, please go to "D. Count Arabidopsis seeds using ImageJ"

# A. Apply our seed counting model in windows

## 1.	To install the tensorflow 1.x, you need the python with version lower than 3.8. 

  * In cmd, type “python” and enter to check if you have the python installed, and what is the version of python. If the version is equal to or higher than 3.8, you need to uninstall the python, and install a lower version. 

  * If your laptop has a python of version lower than 3.8, then you don’t have to do anything.

  * To uninstall your python, do it the same way as you uninstall your windows software.

## 2.	Download the Python 3.7.6 from [here](https://www.python.org/downloads/windows/) by clicking “Windows x86-64 executable installer”

  * When install, select the option “Add Python 3.7 to Path”

## 3.	Install the software needed.  

  * In cmd, type:
  
	`pip install anaconda`
	
	`pip install tensorflow==1.13.2`
	
	`python -m pip install --upgrade pip`
	
	`pip install pillow`
	
	`pip install matplotlib`
	
## 4.	Create a workdir:

  * An example: 
	
	`D:\Projects\Project_done\2022_Seed_count\Test_for_Jupyter`

## 5.	Download the files needed to your work dir

  * Download the file “Files_needed.zip” from [here](https://github.com/ShiuLab/Manuscript_Code/tree/master/2022_Arabidopsis_seed_and_fruit_count) to your workdir, unzip it, and make sure all the folders and files in “Files_needed.zip” are now in your workdir.

## 6.	Count seed using command lines. If you want to use jupyter, please go to Step 7

  * If your scan image contains multiple plate lids, please split the scanned image into images with single plate. Put your images in a folder, here the example is "Scan_images". Run the command below. Don't forget to replace the "D:/Projects/Project_done/2022_Seed_count/Test_for_Jupyter/" with your workdir
	
	`python D:/Projects/Project_done/2022_Seed_count/Test_for_Jupyter/Scripts_for_Faster_R-CNN/00_1_split_scan_images.py D:/Projects/Project_done/2022_Seed_count/Test_for_Jupyter/Scan_images`

  * Go to the image folder, you will find a new folder named "single_plate".

  * Create a folder named “test_images” in your workdir, and put all the single plate seed images you want to count in this folder. Run the code below. Don't forget to replace the "D:/Projects/Project_done/2022_Seed_count/Test_for_Jupyter/" with your workdir. When you see info like "Done. 2022-02-24_12-06-19", the job is done for one image. Go to the folder “test_images” to check the output. 
 
	`python D:/Projects/Project_done/2022_Seed_count/Test_for_Jupyter/Scripts_for_Faster_R-CNN/06_detect_save_image_results.py -base_path D:/Projects/Project_done/2022_Seed_count/Test_for_Jupyter -test_images=test_images`


## 7	Count seed using jupyter

  * Install jupyter

	`python -m pip install jupyter`
	
  * Generate jupyter config file, type:

	`jupyter-notebook --generate-config`

  * Get the path after “Writing default config to:” where the config file is located


  * Change home directory to work directory in config file:

  * Find the config file via the path obtained by Step 4, open it and search for "NotebookApp.notebook_dir", change the path to your workdir

	`NotebookApp.notebook_dir = 'your_work_dir'`


  * Open the jupyter codes by typing:

	`jupyter notebook --notebook-dir=your_workdir\Scripts_for_Faster_R-CNN`

  * example:

	`jupyter notebook --notebook-dir=D:\Projects\Project_done\2022_Seed_count\Test_for_Jupyter\Scripts_for_Faster_R-CNN`


  * Put your images in a folder

  * If your scan image contains multiple plate lids, please split the scanned image into images with single plate. Put your images in a folder, then open the jupyter code "00_1_split_scan_images.py.ipynb". Change the path to your path of this folder with your images. One example is as below. Note that the path should be delimited by “/”, rather than “\”.

	`path = "D:/Projects/Project_done/2022_Seed_count/Test_for_Jupyter/Scan_images"`

  * Run the code "00_1_split_scan_images.py.ipynb". Go to the image folder, you will find a new folder named "single_plate".

  * Create a folder named “test_images” in your workdir, and put all the single plate seed images you want to count in this folder.

  * Open the file "detect_save_image_results.ipynb", and then change the BASE_PATH to your workdir in the file. 

	`BASE_PATH = 'D:/Projects/Project_done/2022_Seed_count/Test_for_Jupyter'`
 
  * Run the code “detect_save_image_results.ipynb”. When you see info like below, the job is done for one image. Go to the folder “test_images” to check the output.
 
	`Done. 2022-02-24_12-06-19`


# B. Apply our seed counting model in Linux

## a. set tensorflow environment

### 1. Install tensorflow(1.X) (version lower than 2.0) in Anaconda 
  * To install TensorFlow (the latest stable release) in a python virtual environment, follow the steps below.
  
  * In Michigan State University HPCC, your may need to load modules as below
  
	`module purge`
	
	`module load GCC/6.4.0-2.28  OpenMPI/2.1.2`
	
	`module load CUDA/10.0.130 cuDNN/7.5.0.56-CUDA-10.0.130`
	
	`module load Python/3.6.4`
	
  * Create a virtual environment
  
	`virtualenv -p python3 tf-1.13.1-env`
	
  * activate this enrironment everytime you use the tensorflow
  
	`source ~/tf-1.13.1-env/bin/activate`
	
  * install tensorflow, for cpu version:
  
	`pip install tensorflow`
	
  * install tensorflow, for gpu version:
  
	`pip install tensorflow-gup==1.13.1`

### 2. Tensorflow object detection API installation

 * Tensorflow version 1.x (version lower than 2.0 is required. Version above 2.0 will not work. If you already have Tensorflow installed and want to check the version of it, please try this: python -c 'import tensorflow as tf; print(tf.__version__)')
* Follow the instruction below to install the API. For the original instruction, check [here](https://github.com/tensorflow/models/blob/master/research/object_detection/g3doc/tf1.md)
* Download models

	`git clone https://github.com/tensorflow/models.git`
	
* Python Package Installation

	cd models/research
	
	`protoc object_detection/protos/*.proto --python_out=.` #NOTE: if can not compile protos, please replace research/object_detection/protos using this compiled [protos](https://github.com/ShiuLab/Manuscript_Code/blob/master/2022_Arabidopsis_seed_and_fruit_count/protos). Or you can download the file "Files_needed.zip" from [here](https://github.com/ShiuLab/Manuscript_Code/tree/master/2022_Arabidopsis_seed_and_fruit_count). Unzip this file, and then you will find the folder "protos" in "research/object_detection". Then using files in this "research/object_detection/protos" to replace files in 
	
	`cp object_detection/packages/tf1/setup.py .`
	
	`python -m pip install .`
  
  * Ensure the Protobuf libraries are compiled and the library directories are added to PYTHONPATH, or set PYTHONPATH in python scripts.
  i.e.:
  
	`import os,sys`
	  
	`sys.path.append("/mnt/home/mengfanr/python-tfgpu/lib/python3.5/site-packages/models")`
	
	`sys.path.append("/mnt/home/mengfanr/python-tfgpu/lib/python3.5/site-packages/models/research")`
	
	`sys.path.append("/mnt/home/mengfanr/python-tfgpu/lib/python3.5/site-packages/models/research/object_detection/utils/")`
	
	`sys.path.append("/mnt/home/mengfanr/python-tfgpu/lib/python3.5/site-packages/models/research/slim/")`

## b. Detect Arabidopsis seeds using our final Faster R-CNN model

### 1. copy files needed
  * Create work directory and copy the files or folders ([graph_train](https://github.com/ShiuLab/Manuscript_Code/tree/master/2022_Arabidopsis_seed_and_fruit_count/Seed_counting_model/graph_train), models/research, [mscoco_label_map.pbtxt](https://github.com/ShiuLab/Manuscript_Code/blob/master/2022_Arabidopsis_seed_and_fruit_count/mscoco_label_map.pbtxt)) to the work directory
  
	`mkdir work_dir`
	
	`cp graph_train models/research mscoco_label_map.pbtxt work_dir -r`

### 2. Seed detection using trained model

  * Create a folder named "test_image", and put your images to be detected within this folder. For testing, images in [Seed_annotation_for_Faster_R-CNN](https://github.com/ShiuLab/Manuscript_Code/tree/master/2022_Arabidopsis_seed_and_fruit_count/Seed_annotation_for_Faster_R-CNN) can be used

	`mkdir test_image`
	
Run the script to count seeds in images

 * `base_path`: the absolute path including the [graph_train](https://github.com/ShiuLab/Manuscript_Code/blob/master/2022_Arabidopsis_seed_and_fruit_count/Seed_counting_model) directory

	`python 06_detect_save_image_results.py --base_path=work_dir --test_images=test_image`
 
 * if you don't want to save image results, you can use:

	`python 06_detect.py --base_path=work_dir --test_images=test_image`


# C. Build your own seed counting models in Linux

  * Here an instruction is given to establish your own seed counting models in Linux

## a. to set the tensorflow environment, do the steps B.a.1 and B.a.2

## b. Training seed detection model

  * The following is for training new/updated Faster-RCNN model.

### 1. Seed annotation

  * This step is only needed for training a new model. For applying the Faster R-CNN model to your seed images, this step is not necessary.
  * During the first round model training, we split one whole plate image into 4 quater images and annotate quater images mannually.

	`python 00_split_images.py image_directory_path`

  * We use [LabelImg](https://github.com/tzutalin/labelImg) (Tzutalin, 2015) to annotate our seeds. LabelImg generates a xml annotation file.
<img src="https://github.com/FanruiMeng/Arabidopsis_seed_count/blob/master/Images/seeds_annotation.png?raw=true"  alt="Seed annotation" height="200"/>
  
### 2. Convert xml files to csv files

  * Script for the conversion is listed below, where annotation is the folder with all the annotation xml files: 

	`pyhton 01_xml_to_csv.py annotation`
  
  * The conversion results in a csv file:

  <table>
  <tr><td><b>filename</b></td> <td><b>width</b></td> <td><b>height</b></td> <td><b>class</b></td> <td><b>x<sub>min<sub></b></td><td><b>y<sub>min<sub></b></td><td><b>x<sub>max<sub></b></td><td><b>y<sub>max<sub></b></td></tr>
  <tr><td>scan21_111918022.jpg</td> <td>2900</td> <td>2900</td> <td>seed</td> <td>1325</td><td>813</td><td>1352</td><td>837</td></tr>
  <tr><td>scan21_111918022.jpg</td> <td>2900</td> <td>2900</td> <td>seed</td> <td>664</td><td>1094</td><td>691</td><td>1116</td></tr>
  <tr><td>..</td> <td>..</td> <td>..</td> <td>..</td> <td>..</td><td>..</td><td>..</td><td>..</td></tr>
  </table>
  
### 3. Convert CSV file to Tensorflow tfrecord file

  * To train the model, the csv files should be converted to tfrecord files.  

	`python 02_generate_tfrecord.py --csv_input=annotation/seeds_labels.csv --output_path=train.record`

### 4. Download tensorflow object detection API pre-trained Faster R-CNN model (inception_v2_coco)

  * In your terminal, issue the following command:

	`wget http://download.tensorflow.org/models/object_detection/faster_rcnn_inception_v2_coco_2018_01_28.tar.gz`

	`tar -xf faster_rcnn_inception_v2_coco_2018_01_28.tar.gz`

### 5. Pipeline configuration

#### Input configuration

  * In the unzip folder `faster_rcnn_inception_v2_coco_2018_01_28`
  * Replace the `pipeline.config` file with [the one provided in this repository](https://github.com/ShiuLab/Manuscript_Code/blob/master/2022_Arabidopsis_seed_and_fruit_count/pipeline.config)
  * For more information on how to change the configuration, see [this document](https://github.com/tensorflow/models/blob/master/research/object_detection/g3doc/configuring_jobs.md).

### 6. Model training

  * Train model with the following parameters (Ren et al., 2017)
    * `logtostderr`: provide logs during training
    * `pipeline_config_path`: the path to the configuration file
    * `train_dir`: output directory for the model
    * `num_clones`: number of processing units used to train the model

	`python 03_train.py --logtostderr --pipeline_config_path=pipeline.config --train_dir=train_dir --num_clones=3`

### 7. Generate a frozen (final) model

  * Generate with the following parameters:
    * `input_type`: the type of input
    * `pipeline_config_path`: the path to the configuration file
    * `trained_checkpoint_prefix`: prefix of the names for the model checkpoint files (ie:"/usr/home/username/train_dir/model.ckpt-10000")
    * `output_directory`: name for output directory

	`python 05_export_inference_graph.py --input_type image_tensor --pipeline_config_path pipeline.config --trained_checkpoint_prefix train_dir/model.ckpt-10000 --output_directory graph_train`

### 8. Detect seeds using trained model

* Parameter: 
  * `base_path`: the absolute path including the [graph_train](https://github.com/ShiuLab/Manuscript_Code/blob/master/2022_Arabidopsis_seed_and_fruit_count/Seed_counting_model) directory
  * `test_images`: test images directory
  
	`python 06_detect.py --base_path=base_path --test_images=test_images`

* if you want to save image results, please use:

	`python 06_detect_save_image_results.py --base_path=base_path --test_images=test_images`

### 9. Accuracy measurement

* Measure accuracy, precision, recall and F1 at IOU 0.5 using 07_01_accuracy_measurement.py <br>

	`python 07_01_accuracy_measurement.py ground.csv detected.csv`

### 10. Estimating seed density

* Determine the average seed number in a circle with a radius of 30 pixels.

	`Rscript seed_density.r`


# D. Count Arabidopsis seeds using ImageJ
   
### 1. Download ImageJ

* Our ImageJ seed count pipeline is developed in Windows.

* The current ImageJ version in windows is 1.53; download it from: [ImageJ](http://wsr.imagej.net/distros/win/ij153-win-java8.zip)  (Schneider et al., 2012) . 

* Make a new folder, `work_dir`; download, save and unzip the ImageJ in `work_dir`.

### 3. Download scripts

* Download these three scripts [here](https://github.com/ShiuLab/Manuscript_Code/tree/master/2022_Arabidopsis_seed_and_fruit_count/ImageJ) to the work_dir: Image_converter.py, small_plate_partial_macro.ijm, seed_image_processing.bat.

### 4. Download test images

* Make a new folder in `work_dir`: `work_dir/images`.

* Some test images can also be found [here](https://github.com/ShiuLab/Manuscript_Code/tree/master/2022_Arabidopsis_seed_and_fruit_count/ImageJ/Example_images), put them in `work_dir/images`.

* For your own images, just put them in the `work_dir/images`.

### 5. Run ImageJ

* Before running, make sure you have installed "Pillow" in your python environment. Otherwise, type "pip install Pillow" in your PC terminal.
* Double click `seed_image_processing.bat`. That's it!

* Note that, when you run `seed_image_processing.bat`, for each image in the "images" folder, the ImageJ interface would be popped up, and when the image is done with seed count, you need to manually close the pop-up ImageJ interface to allow the seed count on next image. Otherwise, it would pause there after one image is done forever.

# References

* Ren SQ, He KM, Girshick R, Sun J. 2017. Faster R-CNN: towards real-time object detection with Region Proposal Networks. Ieee Transactions on Pattern Analysis and Machine Intelligence:1137-1149. https://doi.org/10.1109/TPAMI.2016.2577031.
* Schneider CA, Rasband WS, Eliceiri KW. 2012. NIH Image to ImageJ: 25 years of image analysis. Nature Methods 9:671-675.
* Tzutalin. 2015. LabelImg. Git code. https://github.com/tzutalin/labelImg. [accessed 12 July 2021].