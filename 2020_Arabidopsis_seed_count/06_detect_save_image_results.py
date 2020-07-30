import sys,time
from pathlib import Path
import tensorflow as tf
flags = tf.app.flags
flags.DEFINE_string('base_path', '', 'Path to the CSV input')
flags.DEFINE_string('test_images', '', 'Path to output TFRecord')
FLAGS = flags.FLAGS
from object_detection.utils import label_map_util
from object_detection.utils import visualization_utils as vis_util
import numpy as np
import os
import six.moves.urllib as urllib
import tarfile
import tensorflow as tf
import zipfile
from collections import defaultdict
from io import StringIO
from PIL import Image
from matplotlib import pyplot as plt
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
os.environ['CUDA_VISIBLE_DEVICES'] = "0"
from tensorflow.python.client import device_lib

config = tf.ConfigProto(log_device_placement=True)

#test images dir
PATH_TO_TEST_IMAGES_DIR = FLAGS.test_images
TEST_IMAGE_PATHS=[]
for root, dirs, files in os.walk(PATH_TO_TEST_IMAGES_DIR):
        for f in files:
            if not f.endswith("csv"):
                TEST_IMAGE_PATHS.append(os.path.join(PATH_TO_TEST_IMAGES_DIR,f))    
PATH_TO_CKPT=FLAGS.base_path+"/graph_train/frozen_inference_graph.pb"
PATH_TO_LABELS = FLAGS.base_path+"/mscoco_label_map.pbtxt"
NUM_CLASSES = 1
label_map = label_map_util.load_labelmap(PATH_TO_LABELS)
categories = label_map_util.convert_label_map_to_categories(label_map, max_num_classes=NUM_CLASSES, use_display_name=False)
category_index = label_map_util.create_category_index(categories)


def load_image_into_numpy_array(image):
  (im_width, im_height) = image.size
  return np.array(image.getdata()).reshape(
      (im_height, im_width, 3)).astype(np.uint8)

detection_graph = tf.Graph()
with detection_graph.as_default():
  od_graph_def = tf.GraphDef()
  with tf.gfile.GFile(PATH_TO_CKPT, 'rb') as fid:
    serialized_graph = fid.read()
    od_graph_def.ParseFromString(serialized_graph)
    tf.import_graph_def(od_graph_def, name='')

def run_inference_for_single_image(image, graph):
  with graph.as_default():
    with tf.Session() as sess:
      # Get handles to input and output tensors
      ops = tf.get_default_graph().get_operations()
      all_tensor_names = {output.name for op in ops for output in op.outputs}
      tensor_dict = {}
      for key in ['num_detections', 'detection_boxes', 'detection_scores','detection_classes', 'detection_masks']:
        tensor_name = key + ':0'
        if tensor_name in all_tensor_names:
          tensor_dict[key] = tf.get_default_graph().get_tensor_by_name(tensor_name)
      if 'detection_masks' in tensor_dict:
        # The following processing is only for single image
        detection_boxes = tf.squeeze(tensor_dict['detection_boxes'], [0])
        detection_masks = tf.squeeze(tensor_dict['detection_masks'], [0])
        # Reframe is required to translate mask from box coordinates to image coordinates and fit the image size.
        real_num_detection = tf.cast(tensor_dict['num_detections'][0], tf.int32)
        detection_boxes = tf.slice(detection_boxes, [0, 0], [real_num_detection, -1])
        detection_masks = tf.slice(detection_masks, [0, 0, 0], [real_num_detection, -1, -1])
        detection_masks_reframed = utils_ops.reframe_box_masks_to_image_masks(
            detection_masks, detection_boxes, image.shape[0], image.shape[1])
        detection_masks_reframed = tf.cast(tf.greater(detection_masks_reframed, 0.5), tf.uint8)
        # Follow the convention by adding back the batch dimension
        tensor_dict['detection_masks'] = tf.expand_dims(detection_masks_reframed, 0)
      image_tensor = tf.get_default_graph().get_tensor_by_name('image_tensor:0')
      output_dict = sess.run(tensor_dict,feed_dict={image_tensor: np.expand_dims(image, 0)})
      output_dict['num_detections'] = int(output_dict['num_detections'][0])
      output_dict['detection_classes'] = output_dict['detection_classes'][0].astype(np.uint8)
      output_dict['detection_boxes'] = output_dict['detection_boxes'][0]
      #print(output_dict['detection_boxes'])
      output_dict['detection_scores'] = output_dict['detection_scores'][0]
      if 'detection_masks' in output_dict:
        output_dict['detection_masks'] = output_dict['detection_masks'][0]
  return output_dict
output = open(PATH_TO_TEST_IMAGES_DIR+"/Seed_count_results_"+time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())+".csv","a")
for image_path in TEST_IMAGE_PATHS:
  image = Image.open(image_path)
  (im_width, im_height) = image.size
  image_np = load_image_into_numpy_array(image)
  image_np_expanded = np.expand_dims(image_np, axis=0)
  output_dict = run_inference_for_single_image(image_np, detection_graph)
  scores = output_dict["detection_scores"]
  count = len(scores[scores>=0.5])
  output.write("%s,%d\n"%(image.filename.split("/")[-1],count))
  vis_util.visualize_boxes_and_labels_on_image_array(
      image_np,
      output_dict['detection_boxes'],
      output_dict['detection_classes'],      
      output_dict['detection_scores'],
      category_index,
      instance_masks=output_dict.get('detection_masks'),
      use_normalized_coordinates=True,
      line_thickness=2,
      max_boxes_to_draw=None,
      min_score_thresh=.5,
      skip_scores=True,
      skip_labels=True)
  plt.figure(figsize=((float(format(float(im_width)/float(200),'.2f')), float(format(float(im_height)/float(200),'.2f')))))
  plt.imshow(image_np)
  plt.gca().xaxis.set_major_locator(plt.NullLocator())
  plt.gca().yaxis.set_major_locator(plt.NullLocator())
  plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0) 
  plt.margins(0,0)
  plt.savefig(image_path+"_"+str(count)+"_count.jpg")
