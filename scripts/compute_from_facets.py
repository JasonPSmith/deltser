def run_pydeltser(facets,num_verts):
    dims_dict = [{tuple(facets[i][j]):(j,[]) for j in range(len(facets[i]))} for i in range(len(facets))]
    for k in tqdm.tqdm(range(len(facets)-1,1,-1)):
        for face,boundary in dims_dict[k].items():
            for s in range(k,-1,-1):
                x = tuple(face[:s]+face[s+1:])
                val = dims_dict[k-1].get(x)
                if val != None:
                    boundary[1].append(val[0])
                else:
                    boundary[1].append(len(dims_dict[k-1]))
                    dims_dict[k-1][x] = (len(dims_dict[k-1]),[])
    print('Format built, now calling pydeltser')
    return deltser([[boundary[1] for face,boundary in dims_dict[k].items()] for k in range(len(dims_dict))])
