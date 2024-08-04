# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <Python.h>

# define PY_SSIZE_T_CLEAN
# define EPSILON 0.001
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

void free_memmory_of_vectors_array(int vectors_counted)
{
    int i;
    for(i = 0; i < vectors_counted; i++)
    {
        free(vector_array[i]);
    }
    free(vector_array);
}

void free_memmory_of_centroids(int vectors_counted, double** centroids)
{
    int i;
    for(i = 0; i < vectors_counted; i++)
    {
        free(centroids[i]);
    }
    free(centroids);
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

int initialize_clusters(double** centroids)
{
    int i;
    clusters = (struct cluster *)calloc(sizeof(struct cluster), K);
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
        centroid = (double*)calloc(sizeof(double), vector_length);
        if (centroid == NULL) 
        {
            printf(ERR_MSG);
            free(clusters[i].sum);
            free_memory_of_clusters(i);
            return 1;
        }
        copy_vector_by_cord(centroids[i], centroid);
        clusters[i].centroid = centroid;
        clusters[i].vect_length = vector_length;
        clusters[i].vect_count = 0;
    }
    return 0;
}


void transform_centroids(PyObject *cent_list, double** final_centroids){
    int i, j;
    size_t vect_len, t_i, t_j;
    Py_ssize_t py_len, pos_i, pos_j;
    vect_len = vector_length;
    // Convert size_t to Py_ssize_t
    py_len = (Py_ssize_t)vect_len;
    PyObject *centroid_i, *item_j;
    for (i = 0; i < K; i++)
    {
        t_i = i;
        pos_i = (Py_ssize_t)t_i;
        centroid_i =  PyList_New(py_len);
        for (j = 0; j < vector_length; j++)
        {
            t_j = j;
            pos_j = (Py_ssize_t)t_j;
            item_j = Py_BuildValue("d", final_centroids[i][j]);
            PyList_SetItem(centroid_i, pos_j, item_j);
        }
        PyList_SetItem(cent_list, pos_i, centroid_i);
    }
}

int create_vector_array(PyObject *lst_data_points)
{
    int i, j;
    size_t t_i, t_j;
    Py_ssize_t pos_i, pos_j;
    double value;
    PyObject *item, *vect;
    double* vector_i;
    vector_array = (double**)calloc(sizeof(double*), vectors_count);
    if (vector_array == NULL) 
    {
        printf(ERR_MSG);
        return 1;
    }
    for (i = 0; i < vectors_count; i++)
    {
        t_i = i;
        pos_i = (Py_ssize_t)t_i;
        vector_i = (double *)calloc(sizeof(double), vector_length);
        if (vector_i == NULL) 
        {
            printf(ERR_MSG);
            free_memmory_of_vectors_array(i);
            return 1;
        }
        vect = PyList_GetItem(lst_data_points, pos_i);
        for (j = 0; j < vector_length; j++) 
        {
            t_j = j;
            pos_j = (Py_ssize_t)t_j;
            item = PyList_GetItem(vect, pos_j);
            value = PyFloat_AsDouble(item);
            vector_i[j] = value;
        }
    }
    return 0;
}

int create_initial_centroids(PyObject *lst_centroids)
{
    int i, j;
    size_t t_i, t_j;
    Py_ssize_t pos_i, pos_j;
    double value;
    PyObject *item, *cent;
    double* centroid_i;
    initial_centroids = (double**)calloc(sizeof(double*), K);
    if (initial_centroids == NULL) 
    {
        printf(ERR_MSG);
        return 1;
    }
    for (i = 0; i < K; i++)
    {
        t_i = i;
        pos_i = (Py_ssize_t)t_i;
        centroid_i = (double *)calloc(sizeof(double), K);
        if (centroid_i == NULL) 
        {
            printf(ERR_MSG); 
            free_memmory_of_centroids(i, initial_centroids);
            free_memmory_of_vectors_array(vectors_count);
            return 1;
        }
        cent = PyList_GetItem(lst_centroids, pos_i);
        for (j = 0; j < vector_length; j++) 
        {
            t_j = j;
            pos_j = (Py_ssize_t)t_j;
            item = PyList_GetItem(cent, pos_j);
            value = PyFloat_AsDouble(item);
            centroid_i[j] = value;
        }
    }
    return 0;
}

void turn_clusters_to_centroids(double** final_centroids)
{
    int i;
    for(i = 0; i < K; i++)
    {
        copy_vector_by_cord(clusters[i].centroid, final_centroids[i]);
    }

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
            double *prev_cluster_centroid = calloc(sizeof(double), vector_length);
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
    failure = initialize_clusters(initial_centroids);
    if (failure)
    {
        free_memmory_of_vectors_array(vectors_count);
        return 1;
    }
    free(initial_centroids);
    failure = k_means();
    if (failure)
    {
        free_memmory_of_vectors_array(vectors_count);
        free_memory_of_clusters(K);
        return 1;
    }
    free_memmory_of_vectors_array(vectors_count);
    final_centroids = (double**)calloc(sizeof(double*), K);
    if (final_centroids == NULL) 
    {
        printf(ERR_MSG);
        return 1;
    }
    turn_clusters_to_centroids(final_centroids);
    free_memory_of_clusters(K);
    return 0;
}


int is_double_integer(double value) 
{
    return value == (int)value;
}

static PyObject* fit(PyObject *self, PyObject *args)
{
    PyObject *lst_centroids, *lst_data_points, *cent_list;
    size_t t_k;
    Py_ssize_t size_k;
    if(!PyArg_ParseTuple(args, "iiidoo", &K, &iter_limit, &vector_length, &epsilon, &lst_centroids, &lst_data_points)) 
    {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }
    vectors_count = PyObject_Length(lst_data_points);
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
    failure = activate_kmeans(lst_data_points, initial_centroids);
    if (failure)
    {
        return Py_BuildValue("i", failure);
    }
    t_k = K;
    size_k = (Py_ssize_t)t_k;
    cent_list =  PyList_New(size_k);
    transform_centroids(cent_list, final_centroids);
    free_memmory_of_centroids(K, final_centroids);
    return cent_list; 
}

static PyMethodDef kmeansMethods[] = {
    {"fit",                   /* the Python method name that will be used */
      (PyCFunction) fit, /* the C-function that implements the Python function and returns static PyObject*  */
      METH_VARARGS,           /* flags indicating parameters
accepted for this function */
      PyDoc_STR("A Method that implements the kmeans algorithm")}, /*  The docstring for the function */
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
