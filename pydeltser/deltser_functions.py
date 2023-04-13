import numpy as np
from pydelt import run_deltser

def deltser(face_list, approx=False, approx_val=100000):

    return run_deltser(face_list, approx, approx_val, False, 'null')


def deltser_file(in_file, approx=False, approx_val=100000):

    return run_deltser(np.array([[[]]]), approx, approx_val, True, in_file)

def deltser_facet(in_file, num_vert, approx=False, approx_val=100000):
    print('Loading Facets')
    facets = []
    facet_file = open(in_file)
    for facet in facet_file:
        split_facet = [int(i) for i in facet.strip().split(' ')]
        while len(facets) < len(split_facet):
            facets.append([])
        facets[len(split_facet)-1].append(split_facet)

    print('Computing Faces')
    dims_dict = [{tuple(facets[i][j]):(j,[]) for j in range(len(facets[i]))} for i in range(len(facets))]
    for k in range(len(facets)-1,1,-1):
        for face,boundary in dims_dict[k].items():
            for s in range(k,-1,-1):
                x = tuple(face[:s]+face[s+1:])
                val = dims_dict[k-1].get(x)
                if val != None:
                    boundary[1].append(val[0])
                else:
                    boundary[1].append(len(dims_dict[k-1]))
                    dims_dict[k-1][x] = (len(dims_dict[k-1]),[])

    faces = [[face+(0,) for face in dims_dict[k].keys()] for k in range(len(dims_dict))]
    faces[0] = [(0,) for i in range(num_vert)]

    print('Computing Homology')
    return run_deltser(faces, approx, approx_val, False, 'null')
