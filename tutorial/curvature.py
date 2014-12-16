'''<b>Curvature</b> Curvature measures on the outlines of objects
<hr>
'''

import numpy as np
import scipy.ndimage
from scipy.interpolate import splrep, splev

import cellprofiler.cpmodule as cpm
import cellprofiler.measurements as cpmeas
import cellprofiler.settings as cps
from cellprofiler.cpmath.cpmorphology import get_outline_pts

C_CURVATURE = "Curvature"
FTR_MEAN_CURVATURE = "Mean"
FTR_MEDIAN_CURVATURE = "Median"
FTR_STD_CURVATURE = "StandardDeviation"
FTR_SKEWNESS_CURVATURE = "Skewness"
FTR_KURTOSIS_CURVATURE = "Kurtosis"
ALL_FEATURES = [FTR_MEAN_CURVATURE, FTR_MEDIAN_CURVATURE, FTR_STD_CURVATURE,
                FTR_SKEWNESS_CURVATURE, FTR_KURTOSIS_CURVATURE]

'''Run Example5a using an algorithm that only supports non-overlapping, non-touching'''
SUPPORT_BASIC = "Basic"
'''Run Example5a using an algorithm that supports overlapping'''
SUPPORT_OVERLAPPING = "Overlapping"
'''Run Example5a using an algorithm that supports overlapping and touching'''
SUPPORT_TOUCHING = "Overlapping and touching"

class MeasureCurvature(cpm.CPModule):
    variable_revision_number = 1
    module_name = "MeasureCurvature"
    category = "Measurement"
    
    def create_settings(self):
        self.objects_name = cps.ObjectNameSubscriber("Objects name", "Nuclei")
        self.smoothing = cps.Float(
            "Smoothing factor", 0, 
            doc="""Smoothing factor applied to the spline representation.
            The suggested smoothing factors in the literature are between
            N + A * sqrt(N) where A is the smoothing factor and should
            be betweeen -2 and 2.
            """)
        
    def settings(self):
        return [self.objects_name, self.smoothing]
    
    def run(self, workspace):
        #
        # Get some things we need from the workspace
        #
        measurements = workspace.measurements
        object_set = workspace.object_set
        #
        # Get the objects
        #
        objects_name = self.objects_name.value
        objects = object_set.get_objects(objects_name)
        #
        # Don't look at edges touching the image border or mask
        #
        if objects.has_parent_image:
            mask = objects.parent_image.mask
        else:
            mask = np.ones(objects.shape, bool)
        mask = scipy.ndimage.binary_erosion(
            mask, structure = np.ones((3,3), bool), border_value=False)
        cmean, cmedian, cstd, cskewness, ckurtosis = [
            np.ones(objects.count) * np.nan for _ in range(5)]
        smoothing = self.smoothing.value
        for labels, indexes in objects.get_labels():
            pts, offsets, counts = get_outline_pts(labels, indexes)
            for i in range(len(indexes)):
                if counts[i] < 3:
                    continue
                x = np.arange(counts[i])
                yi, yj = pts[offsets[i]:offsets[i]+counts[i]].transpose()
                s = counts[i] + smoothing * np.sqrt(counts[i])
                si, sj = [splrep(x, y, s=s, per=True) for y in yi, yj]
                (ii, jj), (di, dj), (d2i, d2j) = [[
                    splev(x, s, order) for s in si, sj] for order in 0, 1, 2]
                k = np.abs(di * d2j - d2i * dj) / np.sqrt(di*di + dj*dj)
                k = k[mask[yi, yj]]
                idx = indexes[i] - 1
                cmean[idx] = m = np.mean(k)
                cmedian[idx] = np.median(k)
                cstd[idx] = std = np.std(k)
                kdiff = (k - m) / counts[i]
                cskewness[idx] = np.sum(kdiff * kdiff * kdiff) / std ** 3
                ckurtosis[idx] = np.sum((kdiff * kdiff) * (kdiff * kdiff)) / \
                    std ** 4
        for ftr, values in ((FTR_MEAN_CURVATURE, cmean),
                            (FTR_MEDIAN_CURVATURE, cmedian),
                            (FTR_STD_CURVATURE, cstd),
                            (FTR_SKEWNESS_CURVATURE, cskewness),
                            (FTR_KURTOSIS_CURVATURE, ckurtosis)):
            measurements[self.objects_name.value,
                         self.get_measurement_name(ftr)] = values
            
    def get_measurement_name(self, ftr):
        return "_".join((C_CURVATURE, ftr))
            
    def get_measurement_columns(self, pipeline):
        return [(self.objects_name.value, 
                 self.get_measurement_name(ftr), cpmeas.COLTYPE_FLOAT)
                for ftr in ALL_FEATURES]

    def get_categories(self, pipeline, object_name):
        if object_name == self.objects_name:
            return [C_CURVATURE]
        return []
    
    def get_measurements(self, pipeline, object_name, category):
        if object_name == self.objects_name and category == C_CURVATURE:
            return ALL_FEATURES
        return []