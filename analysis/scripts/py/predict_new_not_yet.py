"""
Melanie Weilert
Stowers Institute
Purpose: Store ancillary functions to predicting using basepairmodels model
"""

import os
import json
import pandas as pd
import numpy as np


def load_basepairmodel(model_dir, use_cpu = True):
    """
    Purpose: Load the keras model.h5 from `basepairmodels`. Requires CustomObject.
    When you use_cpu, then you allocate far less memory into the GPU for predictions.
    """
    import tensorflow as tf
    from tensorflow.keras.utils import CustomObjectScope
    from tensorflow.keras.models import load_model
    from genomicsdlarchsandlosses.bpnet.attribution_prior \
        import AttributionPriorModel
    from genomicsdlarchsandlosses.bpnet.custommodel \
        import CustomModel
    from genomicsdlarchsandlosses.bpnet.losses import \
        MultichannelMultinomialNLL, multinomial_nll, CustomMeanSquaredError

    #Determine whether to use CPU or GPU to compute.
    #Determine whether to use CPU or GPU to compute.
    if use_cpu:
        with tf.device('/cpu:0'):
            with CustomObjectScope({'MultichannelMultinomialNLL': MultichannelMultinomialNLL, 'tf': tf,
                             'CustomMeanSquaredError': CustomMeanSquaredError,
                             'AttributionPriorModel': AttributionPriorModel,
                             'CustomModel': CustomModel}):
#                 os.environ['CUDA_VISIBLE_DEVICES'] = '-1'
                model = load_model(model_dir)
    else:
        with tf.device('/gpu:0'):
            with CustomObjectScope({'MultichannelMultinomialNLL': MultichannelMultinomialNLL, 'tf': tf,
                             'CustomMeanSquaredError': CustomMeanSquaredError,
                             'AttributionPriorModel': AttributionPriorModel,
                             'CustomModel': CustomModel}):
                model = load_model(model_dir)


    return(model)

def predict_basepairmodel(model, seqs, output_width = 1000,
            control_profile = None,
            control_logcount = None,
            return_profile_as_softmax = True,
            return_counts_as_exponentiated = True,
            use_cpu = True,
            batch_size = 128):
    """
    Kudos to Charles: /n/projects/cm2363/bpnet-nucleosomes/work/localimportance/allLocalImportances.py
    Purpose: Given a sequence array, use model to predict outcome.
        + default control_profiles and control_logcounts is to be zeroed out
        + naively detects whether the sequence is one-hot-encoded.
    Input: model object and either a list/array of sequences or a one-hot encoded [regions x position x 4] array.
    Output: [profile array [region x position x tasks], counts array [region x tasks]] per task
    """
    import tensorflow as tf

    #If the sequence is not already one-hot-encoded, then assign./
    seqs = np.array(seqs)
    if (seqs.shape[-1] != 4) and len(seqs.shape)==1:
        seqs = one_hot_encode_sequences(seqs)

    #Create zeroed controls if none are specified
    if not control_profile:
        control_profile = np.zeros((seqs.shape[0], output_width, 2))
    if not control_logcount:
        control_logcount = np.zeros(seqs.shape[0])

    #Define input dictionary for predictions
    input_dict = {'sequence': seqs,
                   'control_profile' : control_profile,
                   'control_logcount' : control_logcount}

    #Determine whether to use CPU or GPU to compute.
    if use_cpu:
        with tf.device('/cpu:0'):
            pred = model.predict(input_dict, batch_size = 1)
    else:
        with tf.device('/gpu:0'):
            #             gpu = tf.config.experimental.list_physical_devices('GPU')
            #             tf.config.experimental.set_memory_growth(gpu[0], True)
            #             tf.config.experimental.set_virtual_device_configuration(gpu[0],
            # [tf.config.experimental.VirtualDeviceConfiguration(memory_limit=4096)])
            import tensorflow_probability as tfp
            pred = model.predict(input_dict, batch_size = batch_size)

    if return_profile_as_softmax:
        pred[0] = pred_logits_to_softmax(pred[0])
    if return_counts_as_exponentiated:
        pred[1] = np.exp(pred[1])
    return(pred)

def one_hot_decode_sequence(array):
    """
    Purpose: Given an array [position x 4], decode sequence to a string.
    """
    onehot_decoder = {
    0: 'A',
    1: 'C',
    2: 'G',
    3: 'T'
    }

    idxs = np.where(array)[1]
    return (''.join([onehot_decoder[i] for i in idxs]))

def one_hot_encode_sequence(sequence):
    """
    Kudos to Charles: /n/projects/cm2363/bpnet-nucleosomes/work/localimportance/allLocalImportances.py
    Purpose: Given a SINGLE sequence string, one-hot encode the data.
        + default control_profiles and control_logcounts is to be zeroed out
        + naively detects whether the sequence is one-hot-encoded.
    """
    onehot_mapping = {
    'A': [1,0,0,0],
    'C': [0,1,0,0],
    'G': [0,0,1,0],
    'T': [0,0,0,1],
    'a': [1,0,0,0],
    'c': [0,1,0,0],
    'g': [0,0,1,0],
    't': [0,0,0,1],
    'N': [0,0,0,0]
    }
    return np.array([onehot_mapping[x] for x in sequence])

def one_hot_encode_sequences(sequences):
    """
    Purpose: Given an array of sequences, one-hot-encode into a [region x position x 4] array.
    """
    return(np.stack([one_hot_encode_sequence(s) for s in sequences]))

def softmax(array):
    """
    Kudos: https://machinelearningmastery.com/softmax-activation-function-with-python/
    """
    assert len(array.shape)==1, 'Input array is not a 1D vector.'
    e = np.exp(array)
    return(e / np.sum(e))

def pred_logits_to_softmax(preds_logits):
    """
    Purpose: Convert [regions x position x channel] logits output of `predict_basepairmodel` to softmax
    """
    assert len(preds_logits.shape)==3, 'Input array is not a 3D tensor [regions x position x task].'
    preds_softmax = np.empty(preds_logits.shape)
    for i in range(preds_logits.shape[0]):
        for j in range(preds_logits.shape[2]):
            preds_softmax[i, :, j] = softmax(preds_logits[i, :, j])
    return(preds_softmax)
