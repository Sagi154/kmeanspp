#define PY_SSIZE_T_CLEAN
#include <Python.h>

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# define ERR_MSG "An Error Has Occurred\n"
# define DBL_MAX 1.7976931348623157e+308

struct cluster
{
    int vect_count;
    int vect_length;
    double* centroid;
    double* sum;
} ;

int vectors_count, vector_length, iter_limit, K, failure;
double epsilon;
double** vector_array;
double** initial_centroids;
double** final_centroids;
struct cluster* clusters = NULL;

void free_memory_of_vectors_array(int vectors_counted)
{
    int i;
    for(i = 0; i < vectors_counted; i++)
    {
        free(vector_array[i]);
    }
    free(vector_array);
}

void free_memory_of_initial_centroids(int vectors_counted)
{
    int i;
    for(i = 0; i < vectors_counted; i++)
    {
        free(initial_centroids[i]);
    }
    free(initial_centroids);
}

void free_memory_of_final_centroids(int vectors_counted)
{
    int i;
    for(i = 0; i < vectors_counted; i++)
    {
        free(final_centroids[i]);
    }
    free(final_centroids);
}

void free_memory_of_clusters(int clusters_count)
{
    int i;
    for(i = 0; i < clusters_count; i++)
    {
        free(clusters[i].centroid);
        free(clusters[i].sum);
    }
    free(clusters);
}

void copy_vector_by_cord(double* copy_from, double* copy_to)
{
    int i;
    for(i = 0; i < vector_length; i++)
    {
        copy_to[i] = copy_from[i];
    }
}

int initialize_clusters()
{
    int i;
    clusters = (struct cluster *)calloc(K, sizeof(struct cluster));
    if (clusters == NULL) 
    {
        printf(ERR_MSG);
        return 1;
    }
    for(i = 0; i < K; i++)
    {
        double* centroid;
        clusters[i].sum = calloc(vector_length, sizeof(double));
        if (clusters[i].sum == NULL) 
        {
            printf(ERR_MSG);
            free_memory_of_clusters(i);
            return 1;
        }
        centroid = (double*)calloc(vector_length, sizeof(double));
        if (centroid == NULL) 
        {
            printf(ERR_MSG);
            free(clusters[i].sum);
            free_memory_of_clusters(i);
            return 1;
        }
        copy_vector_by_cord(initial_centroids[i], centroid);
        clusters[i].centroid = centroid;
        clusters[i].vect_length = vector_length;
        clusters[i].vect_count = 0;
    }
    return 0;
}


int transform_centroids(PyObject *cent_list, double** final_centroids){
    int i, j;
    PyObject* centroid_i;
    PyObject* item_j;
    for (i = 0; i < K; i++)
    {
        centroid_i =  PyList_New(vector_length);
        if (!PyList_Check(centroid_i)) 
        {
            return 1;
        }
        for (j = 0; j < vector_length; j++)
        {
            item_j = Py_BuildValue("d", final_centroids[i][j]);
            if (!item_j) 
            {
                return 1;
            }
            PyList_SetItem(centroid_i, j, item_j);
        }
        PyList_SetItem(cent_list, i, centroid_i);
    }
    return 0;
}

int create_vector_array(PyObject *lst_data_points)
{
    int i, j;
    double value;
    PyObject* item;
    PyObject* vect;
    double* vector_i;
    vector_array = (double**)calloc(vectors_count, sizeof(double*));
    if (vector_array == NULL) 
    {
        printf(ERR_MSG);
        return 1;
    }
    for (i = 0; i < vectors_count; i++)
    {
        vector_i = (double *)calloc(vector_length, sizeof(double));
        if (vector_i == NULL) 
        {
            printf(ERR_MSG);
            free_memory_of_vectors_array(i);
            return 1;
        }
        vect = PyList_GetItem(lst_data_points, i);
        for (j = 0; j < vector_length; j++) 
        {
            item = PyList_GetItem(vect, j);
            value = PyFloat_AsDouble(item);
            vector_i[j] = value;
        }
        vector_array[i] = vector_i;
    }
    return 0;
}

int create_initial_centroids(PyObject *lst_centroids)
{
    int i, j;
    double value;
    PyObject *item;
    PyObject *cent;
    initial_centroids = (double**)calloc(K, sizeof(double*));
    if (initial_centroids == NULL) 
    {
        printf(ERR_MSG);
        return 1;
    }
    for (i = 0; i < K; i++)
    {
        initial_centroids[i] = (double *)calloc(vector_length, sizeof(double));
        if (initial_centroids[i] == NULL) 
        {
            printf(ERR_MSG); 
            free_memory_of_initial_centroids(i);
            free_memory_of_vectors_array(vectors_count);
            return 1;
        }
        cent = PyList_GetItem(lst_centroids, i);
        if (!PyList_Check(cent)) 
        {
            printf(ERR_MSG); 
            free_memory_of_initial_centroids(i);
            free_memory_of_vectors_array(vectors_count);
            return 1;
        }
        for (j = 0; j < vector_length; j++) 
        {
            item = PyList_GetItem(cent, j);
            if (!PyFloat_Check(item)) 
            {
                printf(ERR_MSG); 
                free_memory_of_initial_centroids(i);
                free_memory_of_vectors_array(vectors_count);
                return 1;
            }
            value = PyFloat_AsDouble(item);
            initial_centroids[i][j] = value;
        }
    }
    return 0;
}

int turn_clusters_to_centroids()
{
    int i;
    for(i = 0; i < K; i++)
    {
        final_centroids[i] = (double *)calloc(vector_length, sizeof(double));
        if (final_centroids[i] == NULL) 
        {
            printf(ERR_MSG); 
            free_memory_of_final_centroids(i);
            return 1;
        }
        copy_vector_by_cord(clusters[i].centroid, final_centroids[i]);
    }
    return 0;
}


double calculate_euclidean_distance(double* first_vector, double* second_vector)
{
    double sum_of_points;
    int i;
    sum_of_points = 0.0;
    for (i = 0 ; i < vector_length; i++)
    {
        sum_of_points += pow(second_vector[i] - first_vector[i], 2);
    }
    return sqrt(sum_of_points);
}

int calculate_closest_cluster(struct cluster *clusters, double* vect_xi)
{
    double min_eucledian_distance;
    int min_cluster_index, i;
    min_eucledian_distance = DBL_MAX;
    min_cluster_index = 0;
    for (i = 0 ; i < K; i++) 
    {
        double temp_ed;
        temp_ed = calculate_euclidean_distance(clusters[i].centroid, vect_xi);
        if (temp_ed < min_eucledian_distance)
        {
            min_eucledian_distance = temp_ed;
            min_cluster_index = i;
        }
    }
    return min_cluster_index;
}

void add_xi_to_cluster(struct cluster *cluster, double* vect_xi)
{   
    int i;
    cluster->vect_count = cluster->vect_count + 1;
    for(i = 0; i < vector_length; i++)
    {
        cluster->sum[i] += vect_xi[i];
    }
}

void update_centroid(struct cluster *cluster)
{
    int i;
    for(i = 0; i < vector_length; i++)
    {
        cluster->centroid[i] = cluster->sum[i] / cluster->vect_count;
    }
}

void reset_sum_and_size(struct cluster *cluster)
{
    int i;
    for(i = 0; i < vector_length; i++)
    {
        cluster->sum[i] = 0.0;
    }
    cluster->vect_count = 0;
}



void copy_centroid_to_prev(double* curr_centroid, double* prev_centroid)
{
    int i;
    for (i = 0; i < vector_length; i++)
    {
        prev_centroid[i] = curr_centroid[i];
    }
}

int k_means()
{
    int min_cluster_index, flag, iter_number;
    flag = 0;
    iter_number = 0;
    while(iter_number <= iter_limit && !flag)
    {   
        int i, j;
        iter_number += 1;
        for (i = 0; i < vectors_count; i++)
        {
            min_cluster_index = calculate_closest_cluster(clusters, vector_array[i]);
            add_xi_to_cluster(&clusters[min_cluster_index], vector_array[i]);
        }
        flag = 1;
        for (j = 0; j < K; j++)
        {
            double *prev_cluster_centroid = calloc(vector_length, sizeof(double));
            if (prev_cluster_centroid == NULL) 
            {
                printf(ERR_MSG);
                return 1;
            }
            copy_vector_by_cord(clusters[j].centroid, prev_cluster_centroid);
            update_centroid(&clusters[j]);
            if (flag)
            {
                double convergence = calculate_euclidean_distance(clusters[j].centroid, prev_cluster_centroid);
                if (convergence >= epsilon)
                {
                    flag = 0;
                }
            }
            free(prev_cluster_centroid);
            reset_sum_and_size(&clusters[j]);
        }
    }
    return 0;
}

int activate_kmeans()
{
    failure = initialize_clusters();
    if (failure)
    {
        free_memory_of_initial_centroids(K);
        free_memory_of_vectors_array(vectors_count);
        return 1;
    }
    free_memory_of_initial_centroids(K);
    failure = k_means();
    if (failure)
    {
        free_memory_of_vectors_array(vectors_count);
        free_memory_of_clusters(K);
        return 1;
    }
    free_memory_of_vectors_array(vectors_count);
    final_centroids = (double**)calloc(K, sizeof(double*));
    if (final_centroids == NULL) 
    {
        printf(ERR_MSG);
        free_memory_of_clusters(K);
        return 1;
    }
    failure = turn_clusters_to_centroids();
    if (failure)
    {
        free_memory_of_clusters(K);
        return 1;
    }
    free_memory_of_clusters(K);
    return 0;
}


int is_double_integer(double value) 
{
    return value == (int)value;
}

void print_centroids()
{
    int i, j;
    for (i = 0; i < K; i++)
    {
        for (j = 0; j < vector_length - 1; j++)
        {
            printf("%.4f,", clusters[i].centroid[j]);
        }
        printf("%.4f\n", clusters[i].centroid[vector_length-1]);
    }
}

static PyObject* fit(PyObject *self, PyObject *args)
{
    /**
 * @brief Fits the k-means clustering algorithm to the given data points and centroids.
 *
 * This function initializes and runs the k-means clustering algorithm using the provided data points and initial centroids.
 * It parses the input arguments, processes the data, and returns the final centroids of the clusters.
 *
 * @param self The reference to the module (not used in this function).
 * @param args A tuple containing the following parameters:
 *             - int K: The number of clusters.
 *             - int iter_limit: The maximum number of iterations for the algorithm.
 *             - int vector_length: The length of each data vector.
 *             - double epsilon: The convergence criterion for the algorithm.
 *             - list lst_centroids: A list of initial centroids, where each centroid is represented as a list of doubles.
 *             - list lst_data_points: A list of data points, where each data point is represented as a list of doubles.
 *
 * @return A list of final centroids, where each centroid is represented as a list of doubles. If an error occurs,
 *         it returns an integer indicating the failure status:
 *         - 0: Success
 *         - 1: Error in processing
 *
 * @note This function expects to receive valid input types and dimensions as described. The input lists should be
 *       structured correctly, and all values should be of the expected types to avoid errors.
 */
    PyObject *lst_centroids;
    PyObject *lst_data_points;
    PyObject *cent_list;
    if(!PyArg_ParseTuple(args, "iiidOO", &K, &iter_limit, &vector_length, &epsilon, &lst_centroids, &lst_data_points)) 
    {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }
    vectors_count = PyList_Size(lst_data_points);
    failure = create_vector_array(lst_data_points);
    if (failure)
    {
        return Py_BuildValue("i", failure);
    }
    failure = create_initial_centroids(lst_centroids);
    if (failure)
    {
        return Py_BuildValue("i", failure);
    }
    failure = activate_kmeans();
    if (failure)
    {
        return Py_BuildValue("i", failure);
    }
    cent_list =  PyList_New(K);
    if (!PyList_Check(cent_list)) 
    {
        free_memory_of_final_centroids(K);
        return Py_BuildValue("i", 1);
    }
    failure = transform_centroids(cent_list, final_centroids);
    if (failure)
    {
        free_memory_of_final_centroids(K);
        return Py_BuildValue("i", failure);
    }
    free_memory_of_final_centroids(K);
    return cent_list; 
}

static PyMethodDef kmeansMethods[] = {
    {
        "fit",                   /* the Python method name that will be used */
        fit, /* the C-function that implements the Python function and returns static PyObject*  */
        METH_VARARGS,           /* flags indicating parameters accepted for this function */
        "Fit the k-means clustering algorithm to the given data points and centroids.\n"
        "\n"
        "This function initializes and runs the k-means clustering algorithm using the provided data points and initial centroids.\n"
        "It parses the input arguments, processes the data, and returns the final centroids of the clusters.\n"
        "\n"
        "@param self The reference to the module (not used in this function).\n"
        "@param args A tuple containing the following parameters:\n"
        "             - int K: The number of clusters.\n"
        "             - int iter_limit: The maximum number of iterations for the algorithm.\n"
        "             - int vector_length: The length of each data vector.\n"
        "             - double epsilon: The convergence criterion for the algorithm.\n"
        "             - list lst_centroids: A list of initial centroids, where each centroid is represented as a list of doubles.\n"
        "             - list lst_data_points: A list of data points, where each data point is represented as a list of doubles.\n"
        "\n"
        "@return A list of final centroids, where each centroid is represented as a list of doubles. If an error occurs,\n"
        "         it returns an integer indicating the failure status:\n"
        "         - 0: Success\n"
        "         - 1: Error in processing\n"
        "\n"
        "@note This function expects to receive valid input types and dimensions as described. The input lists should be\n"
        "      structured correctly, and all values should be of the expected types to avoid errors.\n"}, /*  The docstring for the function */
        {NULL, NULL, 0, NULL}     /* The last entry must be all NULL as shown to act as a
                                 sentinel. Python looks for this entry to know that all
                                 of the functions for the module have been defined. */
};

static struct PyModuleDef kmeansmodule = {
    PyModuleDef_HEAD_INIT,
    "mykmeanssp", /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,  /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    kmeansMethods /* the PyMethodDef array from before containing the methods of the extension */
};

PyMODINIT_FUNC PyInit_mykmeanssp(void)
{
    PyObject *m;
    m = PyModule_Create(&kmeansmodule);
    if (!m) {
        return NULL;
    }
    return m;
}

int main(int argc, char **argv) 
{
    return 0;
}
