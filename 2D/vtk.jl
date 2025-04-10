function to_vtk(fname, u)
    interp = vandermonde(rd.element_type, rd.N, equi_nodes(rd.element_type, rd.N)...) / rd.VDM
    pdata = [interp * getindex.(u, 1), interp * pressure.(u, equations)]
    num_elems = ceil(Int, sqrt(md.num_elements ÷ 2))
    vtu_name = MeshData_to_vtk(md, rd, pdata, ["rho", "p"], fname, true)
end

function to_vtk(fname)
    for i in 1:length(sol.t)
        to_vtk(fname * "/step_$i", sol.u[i])
    end
end