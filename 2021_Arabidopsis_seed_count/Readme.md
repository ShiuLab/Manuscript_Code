# We use two methods to count Arabidopsis seeds: Faster R-CNN and ImageJ

# A. Count Arabidopsis seeds using Tensorflow Faster-RCNN model

If you are interested in using Faster-RCNN to detect Arabidopsis seeds. Please see the following Github page:

https://github.com/FanruiMeng/Arabidopsis-Seed-Detection

## Training seed detection model

The following is for training new/updated Faster-RCNN model.

### 1. Tensorflow object detection API installation

* Tensorflow version 1.x (version 2 will not work)
  * [Installation instruction](https://github.com/tensorflow/models/blob/master/research/object_detection/g3doc/tf1.md)
  * Ensure the Protobuf libraries are compiled and the library directories are added to PYTHONPATH, or set PYTHONPATH in python scripts.
  i.e.:
  
   `sys.path.append("/mnt/home/user/python-tfgpu-1.13/lib/python3.6/site-packages/models/research")`
   `sys.path.append("/mnt/home/user/python-tfgpu-1.13/lib/python3.6/site-packages/models/research/object_detection/utils")`
   `sys.path.append("/mnt/home/user/python-tfgpu-1.13/lib/python3.6/site-packages/models/research/slim")`
### 2. Seed annotation

* This is only for training a new model. For applying the Faster R-CNN model, this is not necessary.
* We use [LabelImg](https://github.com/tzutalin/labelImg) to annotate our seeds. LabelImg generate a xml annotation file

<img src="https://github.com/FanruiMeng/Arabidopsis_seed_count/blob/master/Images/seeds_annotation.png?raw=true"  alt="Seed annotation" height="200"/>

* During the first round model training, we split one whole plate image into 4 quater images and annotate quater images mannually.
  
`python 00_split_images.py image_directory_path`
  
### 3. Xml file transform to csv file

* Conversion script: 

`pyhton 01_xml_to_csv.py annotation`
  
* The conversion results in a csv file:

  <table>
  <tr><td><b>filename</b></td> <td><b>width</b></td> <td><b>height</b></td> <td><b>class</b></td> <td><b>xmin</b></td><td><b>ymin</b></td><td><b>xmax</b></td><td><b>ymax</b></td></tr>
  <tr><td>scan21_111918022.jpg</td> <td>2900</td> <td>2900</td> <td>seed</td> <td>1325</td><td>813</td><td>1352</td><td>837</td></tr>
  <tr><td>scan21_111918022.jpg</td> <td>2900</td> <td>2900</td> <td>seed</td> <td>664</td><td>1094</td><td>691</td><td>1116</td></tr>
  <tr><td>..</td> <td>..</td> <td>..</td> <td>..</td> <td>..</td><td>..</td><td>..</td><td>..</td></tr>
  </table>
  
### 4. Convert CSV file to Tensorflow tfrecord file

* Conversion script: 

`python 02_generate_tfrecord.py --csv_input=annotation/seeds_labels.csv --output_path=train.record`

### 5. Download tensorflow object detection API pre-trained Faster R-CNN model

* In your terminal, issue the following command:

`wget http://download.tensorflow.org/models/object_detection/faster_rcnn_inception_v2_coco_2018_01_28.tar.gz`

`tar -xf faster_rcnn_inception_v2_coco_2018_01_28.tar.gz`

### 6. Pipeline configuration

#### Input configuration

* In the unzip folder `faster_rcnn_inception_v2_coco_2018_01_28`
* Replace the `pipeline.config` file with [the one provided in this repository](https://github.com/ShiuLab/Manuscript_Code/blob/master/2020_Arabidopsis_seed_count/pipeline.config)
* For more information on how to change the configuration, see [this document](https://github.com/tensorflow/models/blob/master/research/object_detection/g3doc/configuring_jobs.md).

### 7. Model training

* Train model with the following parameters
  * `logtostderr`: provide logs during training
  * `pipeline_config_path`: the path to the configuration file
  * `train_dir`: output directory for the model
  * `num_clones`: number of processing units used to train the model

`python 03_train.py --logtostderr --pipeline_config_path=pipeline.config --train_dir=train_dir --num_clones=3`

### 8. Generate a frozen (final) model

* Generate with the following parameters:
  * `input_type`: the type of input
  * `pipeline_config_path`: the path to the configuration file
  * `trained_checkpoint_prefix`: prefix of the names for the model checkpoint files (ie:"/usr/home/username/train_dir/model.ckpt-10000")
  * `output_directory`: name for output directory

`python 05_export_inference_graph.py --input_type image_tensor --pipeline_config_path pipeline.config --trained_checkpoint_prefix train_dir/model.ckpt-10000 --output_directory graph_train`

### 9. Detect seeds using trained model

* Parameter: 
  * `base_path`: the absolute path including graph_train directory
  * `test_images`: test images directory
  
`python 06_detect.py --base_path=base_path --test_images=test_images`
* if you want to save image results, please use:

`python 06_detect_save_image_results.py --base_path=base_path --test_images=test_images`

### 10. Accuracy measurement

* Measure accuracy, precision, recall and f1 at IOU 0.5 using 07_01_accuracy_measurement.py <br>

`python 07_01_accuracy_measurement.py ground.csv detected.csv`

### 11. Estimating seed density

* Determine the average seed number in a circle with a radius of 30 pixels.

`Rscript seed_density.r`

# B. Count Arabidopsis seeds using ImageJ
   
### 1. Download ImageJ
* Our ImageJ seed count pipeline is developed in Windows.

* The current ImageJ version in windows is 1.53; download it from: [ImageJ](http://wsr.imagej.net/distros/win/ij153-win-java8.zip). 

* Make a new folder, `work_dir`; download, save and unzip the ImageJ in `work_dir`.

### 3. Download scripts
* Download these three scripts [here](https://github.com/ShiuLab/Manuscript_Code/tree/master/2020_Arabidopsis_seed_count/ImageJ) to the work_dir: Image_converter.py, small_plate_partial_macro.ijm, seed_image_processing.bat.

### 4. Download test images
* Make a new folder in `work_dir`: `work_dir/images`.

* Some test images can also be found [here](https://github.com/ShiuLab/Manuscript_Code/tree/master/2020_Arabidopsis_seed_count/ImageJ/images), put them in `work_dir/images`.

* For your own images, just put them in the `work_dir/images`.

### 5. Run ImageJ
* Before running, make sure you have installed "Pillow" in your python environment. Otherwise, type "pip install Pillow" in your PC terminal.
* Double click `seed_image_processing.bat`. That's it!

* Note that, when you run `seed_image_processing.bat`, for each image in the "images" folder, the ImageJ interface would be popped up, and when the image is done with seed count, you need to manually close the pop-up ImageJ interface to allow the seed count on next image. Otherwise, it would pause there after one image is done forever.
