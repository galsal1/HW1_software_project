import math
import sys

def Kmeans(K, N, d, iter, raw_data):
    assert int(K) > 1 and int(K) < int(N) and K.isdigit(), "Invalid number of clusters!"
    assert int(N) > 1 and N.isdigit(), "Invalid number of points!"
    assert d.isdigit(), "Invalid dimension number of point!"
    assert int(iter) > 1 and int(iter) < 1000 and iter.isdigit(), "Invalid maximum iteration"
    K, N, d, iter = int(K), int(N), int(d), int(iter)
    data = parse_data(raw_data,d,N)
    K_centroids,K_Groups,places = init_classification(K,data,d)
    for cnt in range(iter):
        move_flag = False
        for i in range(len(data)):
            move_flag = move_to_closest(data[i], K_centroids, K_Groups, i, places) or move_flag
        deltasum = 0
        for i in range(K):
            res = groupmean(K_Groups[i], d)
            deltasum += euclidean_distance(res, K_centroids[i])
            K_centroids[i] = res
        if (deltasum < 0.001 or not move_flag):
            break
    return K_centroids

def euclidean_distance(vector1, vector2):
    square_sum = 0
    for i in range(len(vector1)):
        square_sum += math.pow(vector1[i] - vector2[i], 2)
    return math.sqrt(square_sum)

def init_classification(K, data, d):
    K_centroids=[]
    K_Groups=[[] for i in range(K)]
    places=[0 for i in range(len(data))]
    for i in range(K):
        K_centroids.append(data[i])
    for i in range(len(data)):
        vector = data[i]
        distances = [euclidean_distance(vector, K_centroids[j]) for j in range(len(K_centroids))]
        min_dist = min(distances)
        min_index = distances.index(min_dist)
        K_Groups[min_index].append(vector)
        places[i] = min_index
    for i in range(K):
        K_centroids[i] = groupmean(K_Groups[i], d)
    return K_centroids,K_Groups,places

def groupmean(group, d):
    meanvector=[]
    for i in range(d):
        sum=0
        for vector in group:
            sum+=vector[i]
        meanvector.append(sum/len(group))
    return meanvector

def move_to_closest(vector, K_centroids, K_Groups, i ,places):
    current_group = places[i]
    distances = [euclidean_distance(vector, K_centroids[i]) for i in range(len(K_centroids))]
    min_dist = min(distances)
    min_index = distances.index(min_dist)
    if min_index == current_group:
        return False
    else:
        K_Groups[current_group].remove(vector)
        K_Groups[min_index].append(vector)
        places[i] = min_index
        return True

def parse_data(raw_data, d, N):
    raw_data = open(raw_data, "r")
    raw_data_lines = [line.rstrip() for line in raw_data]
    data = []
    for line in raw_data_lines:
        data.append([float(x) for x in line.split(",")])
    return data

if __name__ == '__main__':
    args = sys.argv[1:]
    if len(args) == 5:
        res = Kmeans(args[0], args[1], args[2], args[3], args[4])
    else:
        res = Kmeans(args[0], args[1], args[2], "200", args[3])
    for i in range(int(args[0])):
        print(", ".join("%.4f" % x for x in res[i]))

    print()

    output = open("output_3.txt", "r")
    lines = [line.rstrip() for line in output]
    for line in lines:
        print(line)