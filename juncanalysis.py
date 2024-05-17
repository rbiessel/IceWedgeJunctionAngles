
import numpy as np
import pandas as pd
import geopandas as gpd


import ast


def extract_angles(gdf):

    def str2lst(s):
        s = s.replace('[', '')
        s = s.replace(']', '')
        s = s.split(' ')
        s = [a for a in s if a != '']
        s = np.array(s).astype(np.float32)
        return s

    angles = gdf['angles'].values
    reconangles = []
    for i in range(len(angles)):
        reconangles.append(str2lst(angles[i]))

    return np.concatenate(reconangles)


def orthog_index(angles):
    '''
        Given an arr ay of angles, compute an index where
        0 refers to mostly 60-30 degree angles
        1 refers to mostly 45-90-180 degree angles
    '''
    
    rad = np.deg2rad(angles)
    return np.cos(2 * rad)**2
    # return np.cos(2 * rad)**2

def calculate_angle(junction, line):
    """ 
    """
    jc = np.array(junction.coords[0])[0:2]
    line = np.array(line.coords).T[0:2]
    reference = jc + [0, 1]

    distance = np.sqrt((jc[0] - line[0])**2 + (jc[1] - line[1])**2)
    nearest = np.argmin(distance)

    if nearest == 0:
        nearest += 1
    
    elif nearest == (len(distance) - 1):
        nearest -= 1

    point = line[:, nearest]
    vec1 = reference - jc
    vec2 = point - jc
    unit_vec1 = vec1 / np.linalg.norm(vec1)
    unit_vec2 = vec2 / np.linalg.norm(vec2)
    deg = np.rad2deg(np.arctan2(np.cross(unit_vec2, unit_vec1), np.dot(unit_vec2, unit_vec1)))

    return deg



def get_junctions(gdf, crs):
    intersections = {}
    lines = gdf.geometry.explode(index_parts=True)
    points = []
    intersections = {}

    for l1 in lines:
        for l2 in lines:
            if l1 != l2 and l1.intersects(l2):
                intersection_point = l1.intersection(l2)
                if intersection_point.geom_type == 'Point':
                    points.append(intersection_point)
                    if intersection_point not in intersections:
                        intersections[intersection_point] = set()
                    intersections[intersection_point].add(l1)
                    intersections[intersection_point].add(l2)

    points = []
    jangles = []
    jn = []

    # Calculate angles for each intersection group
    for point, lines in intersections.items():
        if len(lines) == 3:
            angles = []
            points.append(point)
            for line in lines:
                angle = calculate_angle(point, line)
                angles.append(angle)
            angles = sorted(angles)
            angles2 = np.roll(angles, 1)
            difference = np.exp(1j * np.deg2rad(angles)) * np.exp(1j * np.deg2rad(angles2)).conj()
            deg = np.rad2deg(np.angle(difference))
            for i in range(len(deg)):
                if deg[i] < 0:
                    deg[i] += 360
            
            # if (np.abs(np.sum(deg) - 360) < 1):
            jangles.append(str(deg))
            jn.append(len(deg))


    data = {'geometry': points, 'angles': jangles, 'junctions': jn}
    df = pd.DataFrame(data)
    gdf_points = gpd.GeoDataFrame(df, geometry='geometry', crs=crs)
    return gdf_points