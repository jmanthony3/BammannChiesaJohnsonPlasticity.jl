using Ferrite

function vertexdofs(dh::DofHandler, vertexid::VertexIndex)

    cellid, lvidx = vertexid
    nfields = length(dh.field_names)
    sdh = dh.subdofhandlers[dh.cell_to_subdofhandler[cellid]]
    local_vertex_dofs = Int[]
    for ifield in 1:length(sdh.field_names)
        offset = Ferrite.field_offset(sdh, ifield)
        # field_dim = Ferrite.getfielddim(sdh, ifield) # deprecated
        field_dim = Ferrite.n_components(sdh, ifield)
        _field_ip = sdh.field_interpolations[ifield]
        if _field_ip isa Ferrite.VectorizedInterpolation
            field_ip = _field_ip.ip
        else
            field_ip = _field_ip
        end
        vert = Ferrite.vertexdof_indices(field_ip)[lvidx]
        
        for vdof in vert, d in 1:field_dim 
            push!(local_vertex_dofs, (vdof-1)*field_dim + d + offset)
        end
    end

    dofs = zeros(Int, ndofs_per_cell(dh, cellid))
    celldofs!(dofs, dh, cellid)

    return dofs[local_vertex_dofs]
end

function nodeid_to_vertexindex(grid::Grid, nodeid::Int)
    for (cellid, cell) in enumerate(grid.cells)
        for (i, nodeid2) in enumerate(cell.nodes)
            if nodeid == nodeid2
                return VertexIndex(cellid,i)
            end    
        end
    end
    error("Node $(nodeid) does not belong to any cell")
end

# grid = generate_grid(Quadrilateral, (10,10))
# dh = DofHandler(grid)
# add!(dh, :u, Lagrange{2,RefCube,2}()^2) #displacement dofs
# add!(dh, :p, Lagrange{2,RefCube,1}()) #pressure dofs
# close!(dh)

# #Node id to add force on
# nodeid = 11

# #convert nodeid to VertexIndex (cellid, local-node-id) since we need the 
# # the cell id to extract the dofs
# vertexid = nodeid_to_vertexindex(grid, nodeid)

# dof = vertexdofs(dh, vertexid)
# udofs = dof[1:3]
# #f[udof] = [0.0, 0.0, 1.0] #add force