import sys
import pandas as pd
import numpy as np
import math
import mykmeanssp as c
np.random.seed(1234)

def combine_inputs(file_name1, file_name2):
    data_points1 = pd.read_csv(file_name1, header=None)
    data_points2 = pd.read_csv(file_name2, header=None)

    key_column = 0

    data_points1.columns = ['key'] + [f'file1_col{i}' for i in range(1, len(data_points1.columns))]
    data_points2.columns = ['key'] + [f'file2_col{i}' for i in range(1, len(data_points2.columns))]
    data_points1['key'] = data_points1['key'].astype(int)
    data_points2['key'] = data_points2['key'].astype(int)
    combined_data_points = pd.merge(data_points1, data_points2, on='key', suffixes=('_left', '_right'))
    combined_data_points = combined_data_points.sort_values(by='key', ascending=True)
    return combined_data_points

def kmeans_plus_initialization(data_points, K):
    centroids = []
    centroids_indices = []
    indices = [i for i in range(len(data_points))]
    new_centroid_index = np.random.randint(0, len(data_points) - 1)
    centroids.append(data_points[new_centroid_index].tolist())
    centroids_indices.append(new_centroid_index)

    for i in range(K-1):
        new_centroid_index = find_new_centroid(data_points, centroids, indices)
        centroids.append(data_points[new_centroid_index].tolist())
        centroids_indices.append(new_centroid_index)
    return centroids, centroids_indices


def find_new_centroid(data_points, centroids, indices):
    vectors_sum = 0
    vectors_d = [ 0 for i in range(len(data_points))]
    for j in range(len(data_points)):
        vectors_d[j] = find_vect_distance(data_points[j], centroids)
        vectors_sum += vectors_d[j]

    vectors_prob = [vectors_d[i]/vectors_sum for i in range(len(vectors_d))]
    return np.random.choice(indices, p = vectors_prob)


def find_vect_distance(vector, centroids):
    min_distance = sys.maxsize
    for centroid in centroids:
        tmp_distance = calculate_euclidean_distance(vector, centroid)
        if tmp_distance < min_distance:
            min_distance = tmp_distance
    return min_distance


def validity_check(K, iter_limit, epsilon, file_name1, file_name2):
    try:
        if not float(K) == float(int(float(K))):
            print("Invalid number of clusters!")
            return 0
    except Exception as e:
        print("Invalid number of clusters!")
        return 0
    K = int(float(K))

    try:
        if not float(iter_limit) == float(int(float(iter_limit))):
            print("Invalid maximum iteration!")
            return 0
    except Exception as e:
        print("Invalid maximum iteration!")
        return 0
    iter_limit = int(float(iter_limit))

    try:
        if not float(epsilon) >= 0.0:
            print("Invalid epsilon!")
            return 0
    except Exception as e:
        print("Invalid epsilon!")
        return 0
    if not 1 < iter_limit < 1000:
        print("Invalid maximum iteration!")
        return 0
    data_points = combine_inputs(file_name1, file_name2)
    vectors_count = len(data_points)
    if not 1 < K < vectors_count:
        print("Invalid number of clusters!")
        return 0
    return data_points


def calculate_euclidean_distance(vector1, vector2):
    sum = 0
    for i in range(len(vector1)):
        sum += math.pow(vector1[i] - vector2[i], 2)
    return math.sqrt(sum)


def parse_arguments():
    if len(sys.argv) == 6:
        K = sys.argv[1]
        iter_limit = sys.argv[2]
        epsilon = sys.argv[3]
        file_name1 = sys.argv[4]
        file_name2 = sys.argv[5]
        return K, iter_limit, epsilon, file_name1, file_name2
    elif len(sys.argv) == 5:
        K = sys.argv[1]
        iter_limit = 300
        epsilon = sys.argv[2]
        file_name1 = sys.argv[3]
        file_name2 = sys.argv[4]
        return K, iter_limit, epsilon, file_name1, file_name2


def main():
    K, iter_limit, epsilon, file_name1, file_name2 = parse_arguments()
    try:
        data_points = validity_check(K, iter_limit, epsilon, file_name1, file_name2)
        if isinstance(data_points, int):
            return
        data_points = data_points.iloc[:, 1:].to_numpy()
        K = int(float(K))
        iter_limit = int(float(iter_limit))
        epsilon = float(epsilon)
        initialization_centroids, indices = kmeans_plus_initialization(data_points, K)
        vector_len = len(initialization_centroids[0])
        data_points = data_points.tolist()
        final_centroids = c.fit(K, iter_limit, vector_len, epsilon, initialization_centroids, data_points)
        if final_centroids == 1:
            print("An Error Has Occurred")
        else:
            print(",".join(str(index) for index in indices))
            for centroid in final_centroids:
                print(",".join(str("%.4f" % element) for element in centroid))
    except Exception as e:
        print("An Error Has Occurred")

if __name__ == "__main__":
    main()

