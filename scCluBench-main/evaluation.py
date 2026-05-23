import numpy as np
from munkres import Munkres
from sklearn.metrics.cluster import normalized_mutual_info_score as nmi_score
from sklearn.metrics import adjusted_rand_score as ari_score
from sklearn import metrics

def cluster_acc(y_true, y_pred):
    y_true = y_true - np.min(y_true)
    l1 = list(set(y_true))
    numclass1 = len(l1)
    l2 = list(set(y_pred))
    numclass2 = len(l2)

    ind = 0
    if numclass1 != numclass2:
        for i in l1:
            if i in l2:
                pass
            else:
                y_pred[ind] = i
                ind += 1

    l2 = list(set(y_pred))
    numclass2 = len(l2)

    if numclass1 != numclass2:
        print('error')
        return

    cost = np.zeros((numclass1, numclass2), dtype=int)
    for i, c1 in enumerate(l1):
        mps = [i1 for i1, e1 in enumerate(y_true) if e1 == c1]
        for j, c2 in enumerate(l2):
            mps_d = [i1 for i1 in mps if y_pred[i1] == c2]
            cost[i][j] = len(mps_d)

    m = Munkres()
    cost = cost.__neg__().tolist()
    indexes = m.compute(cost)

    new_predict = np.zeros(len(y_pred))
    for i, c in enumerate(l1):
        c2 = l2[indexes[i][1]]
        ai = [ind for ind, elm in enumerate(y_pred) if elm == c2]
        new_predict[ai] = c

    acc = metrics.accuracy_score(y_true, new_predict)
    # y_true：Like 1d array or label indicator array/sparse matrix (correct) label
    # y_pred：Like a one-dimensional array or label indicator array/sparse matrix predicted labels, returned by the classifier
    f1_macro = metrics.f1_score(y_true, new_predict, average='macro')
    return acc, f1_macro


def eva(y_true, y_pred, epoch=0):
    acc, f1 = cluster_acc(y_true, y_pred)
    nmi = nmi_score(y_true, y_pred)
    ari = ari_score(y_true, y_pred)
    print(epoch, ':acc {:.4f}'.format(acc), ', nmi {:.4f}'.format(nmi), ', ari {:.4f}'.format(ari),
          ', f1 {:.4f}'.format(f1))

    # log = open('/home/Ganyanglan/HXY/zinb_sdnc_local/MAGIC_DSC/Mouse/mouse.txt', 'a')
    # if epoch > 298:
    #     log.write('%d\t%.4f\t%.4f\t%.4f\t%.4f\n' %
    #               (epoch,
    #                acc,
    #                nmi,
    #                ari,
    #                f1)
    #               )


# new
from scipy.optimize import linear_sum_assignment as linear_assignment
from sklearn.metrics import f1_score
from sklearn.metrics.cluster import normalized_mutual_info_score as nmi_score
from sklearn.metrics import adjusted_rand_score as ari_score
from sklearn.metrics import fowlkes_mallows_score, v_measure_score, silhouette_score, accuracy_score
from sklearn.metrics.cluster import homogeneity_score, completeness_score
def best_map(y_true, y_pred):
    """
    https://github.com/jundongl/scikit-feature/blob/master/skfeature/utility/unsupervised_evaluation.py
    Permute labels of y_pred to match y_true as much as possible
    """
    if len(y_true) != len(y_pred):
        print("y_true.shape must == y_pred.shape")
        exit(0)

    label_set = np.unique(y_true)
    num_class = len(label_set)

    G = np.zeros((num_class, num_class))
    for i in range(0, num_class):
        for j in range(0, num_class):
            s = y_true == label_set[i]
            t = y_pred == label_set[j]
            G[i, j] = np.count_nonzero(s & t)

    A = linear_assignment(-G)
    new_y_pred = np.zeros(y_pred.shape)
    for i in range(0, num_class):
        new_y_pred[y_pred == label_set[A[1][i]]] = label_set[A[0][i]]
    return new_y_pred.astype(int), label_set[A[1]], label_set[A[0]]

def evaluation(y_true, y_pred):
    y_pred_, label_original, label_truth = best_map(y_true, y_pred)
    acc = accuracy_score(y_true, y_pred_)
    f1_macro = f1_score(y_true, y_pred_, average='macro')
    # f1_micro = f1_score(y_true, best_map(y_true, y_pred), average='micro')
    nmi = nmi_score(y_true, y_pred, average_method='arithmetic')
    ari = ari_score(y_true, y_pred)
    fmi = fowlkes_mallows_score(y_true, y_pred_)
    v_measure = v_measure_score(y_true, y_pred_)
    hom = homogeneity_score(y_true, y_pred_)
    com = completeness_score(y_true, y_pred_)
    # silhouette = silhouette_score(adata.obsm['X_Embeded_z0.6'], aligned_pred_labels)
    # print('origi label', label_original)
    # print('truth label', label_truth)
    # print('recall', recall_score(y_true, y_pred_, average=None))
    # print('precision', precision_score(y_true, y_pred_, average=None))
    return acc, nmi, ari, f1_macro, fmi, v_measure, hom, com, y_pred_
