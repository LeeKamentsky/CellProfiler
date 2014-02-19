'''test_unmixcolors - test the unmixcolors module

CellProfiler is distributed under the GNU General Public License.
See the accompanying file LICENSE for details.

Copyright (c) 2003-2009 Massachusetts Institute of Technology
Copyright (c) 2009-2014 Broad Institute
All rights reserved.

Please see the AUTHORS file for credits.

Website: http://www.cellprofiler.org
'''


import numpy as np
import unittest
from StringIO import StringIO

import cellprofiler.cpimage as cpi
import cellprofiler.cpmodule as cpm
import cellprofiler.measurements as cpmeas
import cellprofiler.pipeline as cpp
import cellprofiler.objects as cpo
import cellprofiler.workspace as cpw
import cellprofiler.modules.unmixcolors as U

INPUT_IMAGE = "inputimage"
def output_image_name(idx):
    return "outputimage%d" % idx

class TestUnmixColors(unittest.TestCase):
    
    def test_01_01_load_v1(self):
        data = r"""CellProfiler Pipeline: http://www.cellprofiler.org
Version:1
SVNRevision:10268

UnmixColors:[module_num:1|svn_version:\'Unknown\'|variable_revision_number:2|show_window:True|notes:\x5B\x5D]
    Stain count:13
    Color image\x3A:Color
    Image name\x3A:Hematoxylin
    Stain:Hematoxylin
    Red absorbance\x3A:0.5
    Green absorbance\x3A:0.5
    Blue absorbance\x3A:0.5
    Image name\x3A:Eosin
    Stain:Eosin
    Red absorbance\x3A:0.5
    Green absorbance\x3A:0.5
    Blue absorbance\x3A:0.5
    Image name\x3A:DAB
    Stain:DAB
    Red absorbance\x3A:0.5
    Green absorbance\x3A:0.5
    Blue absorbance\x3A:0.5
    Image name\x3A:FastRed
    Stain:Fast red
    Red absorbance\x3A:0.5
    Green absorbance\x3A:0.5
    Blue absorbance\x3A:0.5
    Image name\x3A:FastBlue
    Stain:Fast blue
    Red absorbance\x3A:0.5
    Green absorbance\x3A:0.5
    Blue absorbance\x3A:0.5
    Image name\x3A:MethylGreen
    Stain:Methyl green
    Red absorbance\x3A:0.5
    Green absorbance\x3A:0.5
    Blue absorbance\x3A:0.5
    Image name\x3A:AEC
    Stain:AEC
    Red absorbance\x3A:0.5
    Green absorbance\x3A:0.5
    Blue absorbance\x3A:0.5
    Image name:AnilineBlue
    Stain:Aniline blue
    Red absorbance\x3A:0.5
    Green absorbance\x3A:0.5
    Blue absorbance\x3A:0.5
    Image name\x3A:Azocarmine
    Stain:Azocarmine
    Red absorbance\x3A:0.5
    Green absorbance\x3A:0.5
    Blue absorbance\x3A:0.5
    Image name\x3A:AlicanBlue
    Stain:Alican blue
    Red absorbance\x3A:0.5
    Green absorbance\x3A:0.5
    Blue absorbance\x3A:0.5
    Image name\x3A:PAS
    Stain:PAS
    Red absorbance\x3A:0.5
    Green absorbance\x3A:0.5
    Blue absorbance\x3A:0.5
    Image name\x3A:HematoxylinAndPAS
    Stain:Hematoxylin and PAS
    Red absorbance\x3A:0.5
    Green absorbance\x3A:0.5
    Blue absorbance\x3A:0.5
    Image name\x3A:RedWine
    Stain:Custom
    Red absorbance\x3A:0.1
    Green absorbance\x3A:0.2
    Blue absorbance\x3A:0.3
"""
        pipeline = cpp.Pipeline()
        def callback(caller, event):
            self.assertFalse(isinstance(event, cpp.LoadExceptionEvent))
        pipeline.add_listener(callback)
        pipeline.load(StringIO(data))
        self.assertEqual(len(pipeline.modules()), 1)
        module = pipeline.modules()[0]
        self.assertTrue(isinstance(module, U.UnmixColors))
        self.assertEqual(module.input_image_name, "Color")
        self.assertEqual(module.method, U.METHOD_CHOOSE_STAINS)
        self.assertEqual(module.stain_count.value, 13)
        self.assertEqual(module.outputs[0].image_name, "Hematoxylin")
        self.assertEqual(module.outputs[-1].image_name, "RedWine")
        for i, stain in enumerate((
            U.CHOICE_HEMATOXYLIN, U.CHOICE_EOSIN, U.CHOICE_DAB,
            U.CHOICE_FAST_RED, U.CHOICE_FAST_BLUE, U.CHOICE_METHYL_GREEN,
            U.CHOICE_AEC, U.CHOICE_ANILINE_BLUE, U.CHOICE_AZOCARMINE,
            U.CHOICE_ALICAN_BLUE, U.CHOICE_PAS)):
            self.assertEqual(module.outputs[i].stain_choice, stain)
        self.assertAlmostEqual(module.outputs[-1].red_absorbance.value, .1)
        self.assertAlmostEqual(module.outputs[-1].green_absorbance.value, .2)
        self.assertAlmostEqual(module.outputs[-1].blue_absorbance.value, .3)
        
    def test_01_03_load_v3(self):
        data = r"""CellProfiler Pipeline: http://www.cellprofiler.org
Version:3
DateRevision:20140219141600
GitHash:8e5019f
ModuleCount:2
HasImagePlaneDetails:False

UnmixColors:[module_num:1|svn_version:\'Unknown\'|variable_revision_number:3|show_window:True|notes:\x5B\x5D|batch_state:array(\x5B\x5D, dtype=uint8)|enabled:True]
    Stain count:1
    Select the input color image:Color
    Unmixing method:Two-color unmixing
    Second output image:Eosin1
    First output image:Hematoxylin1
    Stain:Hematoxylin
    Red absorbance:0.4
    Green absorbance:0.5
    Blue absorbance:0.6

UnmixColors:[module_num:2|svn_version:\'Unknown\'|variable_revision_number:3|show_window:True|notes:\x5B\x5D|batch_state:array(\x5B\x5D, dtype=uint8)|enabled:True]
    Stain count:2
    Select the input color image:Color5
    Unmixing method:Choose stains
    Second output image:Eosin
    Name the output name:Hematoxylin1
    Stain:Hematoxylin
    Red absorbance:0.7
    Green absorbance:0.8
    Blue absorbance:0.9
    Name the output name:Eosin1
    Stain:Eosin
    Red absorbance:0.1
    Green absorbance:0.2
    Blue absorbance:0.3

"""
        pipeline = cpp.Pipeline()
        def callback(caller, event):
            self.assertFalse(isinstance(event, cpp.LoadExceptionEvent))
        pipeline.add_listener(callback)
        pipeline.load(StringIO(data))
        self.assertEqual(len(pipeline.modules()), 2)
        module = pipeline.modules()[0]
        self.assertTrue(isinstance(module, U.UnmixColors))
        self.assertEqual(module.input_image_name, "Color")
        self.assertEqual(module.method, U.METHOD_TWO_COLOR_UNMIXING)
        self.assertEqual(module.outputs[0].image_name, "Hematoxylin1")
        self.assertEqual(module.second_output_image, "Eosin1")
        self.assertEqual(module.outputs[0].stain_choice, U.CHOICE_HEMATOXYLIN)
        self.assertEqual(module.outputs[0].red_absorbance, .4)
        self.assertEqual(module.outputs[0].green_absorbance, .5)
        self.assertEqual(module.outputs[0].blue_absorbance, .6)
        
        module = pipeline.modules()[1]
        self.assertTrue(isinstance(module, U.UnmixColors))
        self.assertEqual(module.input_image_name, "Color5")
        self.assertEqual(module.method, U.METHOD_CHOOSE_STAINS)
        self.assertEqual(len(module.outputs), 2)
        for output, image_name, stain_choice, r, g, b in (
            (module.outputs[0], "Hematoxylin1", U.CHOICE_HEMATOXYLIN, .7, .8, .9),
            (module.outputs[1], "Eosin1", U.CHOICE_EOSIN, .1, .2, .3)):
            self.assertEqual(output.image_name, image_name)
            self.assertEqual(output.stain_choice, stain_choice)
            self.assertEqual(output.red_absorbance, r)
            self.assertEqual(output.green_absorbance, g)
            self.assertEqual(output.blue_absorbance, b)
             
    def make_workspace(self, pixels, choices):
        '''Make a workspace for running UnmixColors
        
        pixels - input image
        choices - a list of choice strings for the images desired
        '''
        pipeline = cpp.Pipeline()
        def callback(caller, event):
            self.assertFalse(isinstance(event, cpp.RunExceptionEvent))
        pipeline.add_listener(callback)
        
        module = U.UnmixColors()
        module.input_image_name.value = INPUT_IMAGE
        module.outputs[0].image_name.value = output_image_name(0)
        module.outputs[0].stain_choice.value = choices[0]
        for i, choice in enumerate(choices[1:]):
            module.add_image()
            module.outputs[i+1].image_name.value = output_image_name(i+1)
            module.outputs[i+1].stain_choice.value = choice
        
        module.module_num = 1
        pipeline.add_module(module)
        
        image_set_list = cpi.ImageSetList()
        image_set = image_set_list.get_image_set(0)
        image = cpi.Image(pixels)
        image_set.add(INPUT_IMAGE, image)
        
        workspace = cpw.Workspace(pipeline, module, image_set, cpo.ObjectSet(),
                                  cpmeas.Measurements(), image_set_list)
        return workspace, module
    
    @staticmethod
    def make_image(expected, absorbances):
        eps = 1.0/256.0/2.0
        absorbance = 1 - expected
        log_absorbance = np.log(absorbance + eps)
        absorbances = np.array(absorbances)
        absorbances = absorbances / np.sqrt(np.sum(absorbances**2))
        log_absorbance = log_absorbance[:,:,np.newaxis] * absorbances[np.newaxis, np.newaxis, :]
        image = np.exp(log_absorbance) - eps
        return image
    
    def test_02_01_zeros(self):
        '''Test on an image of all zeros'''
        workspace, module = self.make_workspace(np.zeros((10,20,3)),
                                                [U.CHOICE_HEMATOXYLIN])
        module.run(workspace)
        image = workspace.image_set.get_image(output_image_name(0))
        #
        # All zeros in brightfield should be all 1 in stain
        #
        np.testing.assert_almost_equal(image.pixel_data, 1, 2)
        
    def test_02_02_ones(self):
        '''Test on an image of all ones'''
        workspace, module = self.make_workspace(np.ones((10,20,3)),
                                                [U.CHOICE_HEMATOXYLIN])
        module.run(workspace)
        image = workspace.image_set.get_image(output_image_name(0))
        #
        # All ones in brightfield should be no stain
        #
        np.testing.assert_almost_equal(image.pixel_data, 0, 2)
        
    def test_02_03_one_stain(self):
        '''Test on a single stain'''
        
        np.random.seed(23)
        expected = np.random.uniform(size=(10,20))
        image = self.make_image(expected, U.ST_HEMATOXYLIN)
        workspace, module = self.make_workspace(image, [U.CHOICE_HEMATOXYLIN])
        module.run(workspace)
        image = workspace.image_set.get_image(output_image_name(0))
        np.testing.assert_almost_equal(image.pixel_data, expected, 2)
        
    def test_02_04_two_stains(self):
        '''Test on two stains mixed together'''
        np.random.seed(24)
        expected_1 = np.random.uniform(size=(10,20)) * .5
        expected_2 = np.random.uniform(size=(10,20)) * .5
        #
        # The absorbances should add in log space and multiply in
        # the image space
        #
        image = self.make_image(expected_1, U.ST_HEMATOXYLIN)
        image *= self.make_image(expected_2, U.ST_EOSIN)
        workspace, module = self.make_workspace(image, [
            U.CHOICE_HEMATOXYLIN, U.CHOICE_EOSIN ])
        module.run(workspace)
        image_1 = workspace.image_set.get_image(output_image_name(0))
        np.testing.assert_almost_equal(image_1.pixel_data, expected_1, 2)
        image_2 = workspace.image_set.get_image(output_image_name(1))
        np.testing.assert_almost_equal(image_2.pixel_data, expected_2, 2)
        
    def test_02_05_custom_stain(self):
        '''Test on a custom value for the stains'''
        np.random.seed(25)
        absorbance = np.random.uniform(size=3)
        expected = np.random.uniform(size=(10,20))
        image = self.make_image(expected, absorbance)
        workspace, module = self.make_workspace(image, [U.CHOICE_CUSTOM])
        ( module.outputs[0].red_absorbance.value,
          module.outputs[0].green_absorbance.value,
          module.outputs[0].blue_absorbance.value ) = absorbance
        module.run(workspace)
        image = workspace.image_set.get_image(output_image_name(0))
        np.testing.assert_almost_equal(image.pixel_data, expected, 2)
        
    def test_02_06_two_color(self):
        r = np.random.RandomState()
        r.seed(26)
        intensities = r.uniform(size=(20, 20, 2))
        color1 = np.array(U.ST_HEMATOXYLIN)
        color2 = np.array(U.ST_EOSIN)
        
        img1, img2 = [ 
            c[np.newaxis, np.newaxis, :] +
            (1 - intensity[:, :, np.newaxis]) * (1 - c[np.newaxis, np.newaxis, :])
            for intensity, c in ((intensities[:, :, 0], color1),
                                 (intensities[:, :, 1], color2))]
        img = np.exp(img1 + img2)
        workspace, module = self.make_workspace(img, [U.CHOICE_HEMATOXYLIN])
        assert isinstance(module, U.UnmixColors)
        module.method.value = U.METHOD_TWO_COLOR_UNMIXING
        module.second_output_image.value = output_image_name(1)
        module.run(workspace)
        img1 = workspace.image_set.get_image(output_image_name(0)).pixel_data
        img2 = workspace.image_set.get_image(output_image_name(1)).pixel_data
        
