#include "cellSchAlgCPU.h"
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <chrono>
#include <omp.h>

// #include <pybind11/pybind11.h>
// #include <pybind11/stl.h>

namespace dpnblist
{
    CellList::CellList(const Box &box, float r_cutoff) : _box(box), _r_cutoff(r_cutoff)
    {
        std::array<float, 3> box_length = box.get_lengths_cpu3();
        // _cell_length[0] = box_length[0] / _r_cutoff;
        // _cell_length[1] = box_length[1] / _r_cutoff;
        // _cell_length[2] = box_length[2] / _r_cutoff;
        for (int i = 0; i < 3; i++)
        {
            _cell_length[i] = box_length[i] / _r_cutoff;
            new_cutoff[i] = box_length[i] / _cell_length[i];
        }
        _ncells = _cell_length[0] * _cell_length[1] * _cell_length[2];
        // A = {_cell_length[0], 0, 0, 0, _cell_length[1], 0, 0, 0, _cell_length[2]};
        // inv_A = A.invert();
        offset = {
            {-1, -1, -1},{-1, -1, 0},{-1, -1, 1},{-1, 0, -1},{-1, 0, 0},{-1, 0, 1},{-1, 1, -1},{-1, 1, 0},{-1, 1, 1},
            {0, -1, -1},{0, -1, 0},{0, -1, 1},{0, 0, -1},{0, 0, 1},{0, 1, -1},{0, 1, 0},{0, 1, 1},
            {1, -1, -1},{1, -1, 0},{1, -1, 1},{1, 0, -1},{1, 0, 0},{1, 0, 1},{1, 1, -1},{1, 1, 0},{1, 1, 1}
        };
    }

    int CellList::get_cell_index(std::array<int, 3> &cell_vector) const
    {
        /* calculate the cell_index according to cell_vector
        * c = cxLcyLcz + cyLcz + cz
        */
        return cell_vector[0] * _cell_length[1] * _cell_length[2] + cell_vector[1] * _cell_length[2] + cell_vector[2];
    }

    std::array<int, 3> CellList::get_cell_vector(int cell_index) const
    {
        /* calculate the cell_vector according to cell_index
            * cx = c / (LcyLcz)
            * cy = (c / Lcz) - cxLcy or (c / Lcz) mod Lcy
            * cz = c mod Lcz
        */
        std::array<int, 3> cell_vector;
        cell_vector[0] = cell_index / (_cell_length[1] * _cell_length[2]);
        // cell_vector[1] = std::fmod((cell_index / _cell_length[2]), _cell_length[1]);
        cell_vector[1] = (cell_index / _cell_length[2]) % _cell_length[1];
        cell_vector[2] = cell_index % _cell_length[2];
        return cell_vector;
    }

    void CellList::build(std::vector<std::array<float, 3>> &xyz)
    {
        // build the cell list
        _natoms = xyz.size();
    }

    void CellList::update(std::vector<std::array<float, 3>> &xyz)
    {
        reset();
        // update the linked cell list
        std::array<int, 3> xyz_cell_vec;

        for (int i = 0; i < _natoms; i++)
        {
            for (int j = 0; j < 3; j++){                            // if the cell index is out of range, move it back
                xyz_cell_vec[j] = xyz[i][j] / new_cutoff[j];
                if (xyz_cell_vec[j] == _cell_length[j])
                    xyz_cell_vec[j] = xyz_cell_vec[j] - 1;
            }
            int cell_index = get_cell_index(xyz_cell_vec);     // calculate the cell index
            _count[cell_index]++;
            _lscl[i] = _head[cell_index];                           // update the lscl
            _head[cell_index] = i;                                  // update the head
        }
    }

    void CellList::reset()
    {
        _head.clear();
        _head.resize(_ncells, -1);
        _lscl.clear();
        _lscl.resize(_natoms, -1);
        _count.clear();
        _count.resize(_ncells, 0);
    }

    int CellList::get_ncells() const
    {
        return _ncells;
    }
    /*
    // get the atoms in the cell and return
    std::vector<int> CellList::get_atoms_in_cell(int cell_index) const
    {
        // while (atom_index != EMPTY)
        int n = _count[cell_index];
        std::vector<int> atoms_in_cell(n, 0);
        int atom_index = _head[cell_index];
        for (int i = 0; i < n; ++i)
        {
            atoms_in_cell[i] = atom_index;
            atom_index = _lscl[atom_index];
        }
        return atoms_in_cell;
    }
    */
    std::vector<int> CellList::get_atoms_in_cell(int cell_index) const
    {
        if (cell_index >= get_ncells())
        {
            throw std::out_of_range("cell index out of range");
        }
        std::vector<int> atoms_in_cell;
        int atom_index = _head[cell_index];
        while (atom_index != -1)
        {
            atoms_in_cell.push_back(atom_index);
            atom_index = _lscl[atom_index];
        }
        return atoms_in_cell;
    }
    std::vector<int> CellList::get_head(){
        return _head;
    };
    std::vector<int> CellList::get_lscl(){
        return _lscl;
    };
    std::vector<int> CellList::get_count(){
        return _count;
    };
    /*
        * get the neighbors of the cell
        * the cell is described by a vector, the cell_vector is the vector of the central cell, the matrix is the offset vector, the (matrix + cell_vector) is the 26 neighbors cell vector.
        * the vector can be express as r = f * A, A is the box matrix, f is the fractional coordinate, can be express as f = inv_A * r
        * if periodic boundary condition is considered, wrapped_f = f - floor(f), wrapped_f is the wrapped fractional coordinate, so t = A * wrapped_f, t is the wrapped vector
        * reference: https://scicomp.stackexchange.com/questions/20165/periodic-boundary-conditions-for-triclinic-box
    */
    /*
    std::vector<int> CellList::get_neighbors(int cell_index) const
    {
        std::vector<int> neighbors = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, cell_index};
        Vec3<int> cell_vector = get_cell_vector(cell_index);
        // define a offset matrix 26*3
        // std::vector<Vec3<int>> offset = {
        //     {-1, -1, -1},{-1, -1, 0},{-1, -1, 1},{-1, 0, -1},{-1, 0, 0},{-1, 0, 1},{-1, 1, -1},{-1, 1, 0},{-1, 1, 1},
        //     {0, -1, -1},{0, -1, 0},{0, -1, 1},{0, 0, -1},{0, 0, 1},{0, 1, -1},{0, 1, 0},{0, 1, 1},
        //     {1, -1, -1},{1, -1, 0},{1, -1, 1},{1, 0, -1},{1, 0, 0},{1, 0, 1},{1, 1, -1},{1, 1, 0},{1, 1, 1} 
        // };
        // Mat3<float> A = {_cell_length[0], 0, 0, 0, _cell_length[1], 0, 0, 0, _cell_length[2]};
        // Mat3<float> inv_A = A.invert();
        Vec3<float> f;
        Vec3<float> plus = {0.00001, 0.00001, 0.00001};
        Vec3<float> wrapped_f;
        Vec3<float> wrap_neb_vector;
        Vec3<int> round_wrap_neb_vector;

        for (int i = 0; i < 26; ++i) {
            Vec3<int> neb_vector = cell_vector + offset[i];
            for (int k = 0; k < 3; ++k) {
                f[k] = neb_vector[k] / _cell_length[k];
                wrapped_f[k] = f[k] - std::floor(f[k] + 0.000001);
                wrap_neb_vector[k] = wrapped_f[k] * _cell_length[k];
                round_wrap_neb_vector[k] = round(wrap_neb_vector[k]);
            }
            neighbors[i] = get_cell_index(round_wrap_neb_vector);
        }
        return neighbors;
    }
    */

    // std::vector<int> CellList::get_neighbors(int cell_index, std::vector<int> &neighbors, std::vector<Vec3<int>> &shifts) const
    void CellList::get_neighbors(int cell_index, std::array<int, 27> &neighbors, std::vector<std::array<int, 3>> &shifts) const
    {
        // std::vector<int> neighbors(27, 0);
        std::array<int, 3> cell_vector = get_cell_vector(cell_index);
        int x_neb, y_neb, z_neb;
        int n = 0;
        for (int x = cell_vector[0] - 1; x <= cell_vector[0] + 1; x++)
        {
            for (int y = cell_vector[1] - 1; y <= cell_vector[1] + 1; y++)
            {
                for (int z = cell_vector[2] - 1; z <= cell_vector[2] + 1; z++)
                {
                    x_neb = x;
                    y_neb = y;
                    z_neb = z;
                    std::array<int, 3> shift = {0, 0, 0};
                    // periodic
                    if (x < 0)
                    {
                        x_neb = x + _cell_length[0];
                        shift[0] = -1;
                    }
                    else if (x >= _cell_length[0])
                    {
                        x_neb = x - _cell_length[0];
                        shift[0] = 1;
                    }

                    if (y < 0)
                    {
                        y_neb = y + _cell_length[1];
                        shift[1] = -1;
                    }
                    else if (y >= _cell_length[1])
                    {
                        y_neb = y - _cell_length[1];
                        shift[1] = 1;
                    }

                    if (z < 0)
                    {
                        z_neb = z + _cell_length[2];
                        shift[2] = -1;
                    }
                    else if (z >= _cell_length[2])
                    {
                        z_neb = z - _cell_length[2];
                        shift[2] = 1;
                    }
                    shifts[n] = shift;
                    std::array<int, 3> neighbor_cell_vector = {x_neb, y_neb, z_neb};
                    int neighbor_cell_index = get_cell_index(neighbor_cell_vector);
                    if (neighbor_cell_index != cell_index)
                    {
                        neighbors[n] = neighbor_cell_index;
                        n++;
                    }
                }
            }
        }
        neighbors[n] = cell_index;
        // return neighbors;
    }


    CellSearchCPU::CellSearchCPU(const Box& box, const std::vector<std::vector<float>>& xyz, float r_cutoff): _box(box), _xyz(xyz), _r_cutoff(r_cutoff), _cell_list(box, r_cutoff)
    {
        //float volume = box.get_lengths_cpu3()[0] * box.get_lengths_cpu3()[1] * box.get_lengths_cpu3()[2];
        //float pre_nnebs = xyz.size() / volume * 4 * 3.1415926 * r_cutoff * r_cutoff * r_cutoff / 3; // M_PI
        //_nnebs = int(pre_nnebs * 1.11 + 61);

        std::array<float,3> cube_size = box.get_lengths_cpu3();
        _nnebs = static_cast<int>(std::ceil(1.5 * (5 * xyz.size() / (cube_size[0] * cube_size[1] * cube_size[2])) * r_cutoff * r_cutoff * r_cutoff + 45));
        if (_nnebs < 100) _nnebs = 100;
    }

    CellSearchCPU::~CellSearchCPU()
    {
        delete[] _neighlist_1d;
    }

    void CellSearchCPU::reset()
    {
        _neighborListArray.clear();
        _neighborListArray.resize(_natoms);
    }

    void CellSearchCPU::get_neighbor_cell_array()
    {
        // #pragma omp parallel for
        for (int i = 0; i < _ncells; i++)
        {
            std::array<int, 27> neighbors = {};
            std::vector<std::array<int, 3>> shifts(27, {0, 0, 0});
            _cell_list.get_neighbors(i, neighbors, shifts);
            _neighbor_cell_Array[i] = neighbors;
            _neighbor_cell_Array_shifts[i] = shifts;
        }
    }

    //void CellSearchCPU::search(std::vector<std::vector<float>> &xyz)
    void CellSearchCPU::search()
    {
        _natoms = _xyz.size();
        //std::cout << _natoms << std::endl;
        std::vector<std::array<float, 3>> vec3xyz;
        vec3xyz.resize(_natoms);
        for (int i = 0; i < _natoms; i++)
        {
            vec3xyz[i][0] = _xyz[i][0];
            vec3xyz[i][1] = _xyz[i][1];
            vec3xyz[i][2] = _xyz[i][2];
        }
        build(vec3xyz);
    }
    void CellSearchCPU::build(std::vector<std::array<float, 3>> &xyz)
    {
        _natoms = xyz.size();
        _ncells = _cell_list.get_ncells();
        //std::cout << _ncells << std::endl;

        _neighbor_cell_Array.resize(_ncells);
        _neighbor_cell_Array_shifts.resize(_ncells);
        get_neighbor_cell_array();
        // time milliseconds
        auto cell_start = std::chrono::high_resolution_clock::now();
        _cell_list.build(xyz);
        auto cell_end = std::chrono::high_resolution_clock::now();
        auto cell_duration = std::chrono::duration_cast<std::chrono::milliseconds>(cell_end - cell_start);
        //std::cout << "build cell list time: " << cell_duration.count() << "ms" << std::endl;
        update(xyz);
    }

    /*
        循环遍历每个cell，再循环每个cell中的原子，再遍历其邻居cell，最后遍历邻居cell中的原子，计算距离，如果距离小于截断半径，将其加入邻居列表
        Loop through each cell, then loop through the atoms in the cell, then loop through the neighboring cells, 
        and finally loop through the atoms in the neighboring cells, calculate the distance, 
        if the distance is less than the cutoff radius, add it to the neighbor list
    */
    void CellSearchCPU::update(std::vector<std::array<float, 3>> &xyz)
    {
        reset();
        //bool is_update = judge_update(xyz);
        bool is_update = true;
        if (is_update){
            _box.wrap_p(xyz);
            _cell_list.update(xyz);
            std::vector<int> _head = _cell_list.get_head();
            std::vector<int> _lscl = _cell_list.get_lscl();
            std::vector<int> _count = _cell_list.get_count();
            
            auto build_start = std::chrono::high_resolution_clock::now();
            _neighlist_1d = new int[_nnebs * _natoms];
            //if (_neighlist_1d == nullptr){
                //_neighlist_1d = new int[_nnebs * _natoms];
            //}
            for (int i = 0; i < _nnebs * _natoms; i++)
            {
                _neighlist_1d[i] = -1;
            }
            float rc2 = _r_cutoff * _r_cutoff;
            std::array<float,3> box_length = _box.get_lengths_cpu3();
            #pragma omp parallel for
            for (int cell_index = 0; cell_index < _ncells; ++cell_index)                 // Loop through each cell
            // for (int cell_index = 0; cell_index < 1; ++cell_index)                 // Loop through each cell
            {
                // int nebcout = 0;
                int ni = _count[cell_index];// 查找单元内有几个原子
                int i = _head[cell_index];// 根据单元序号查找单元的第一个原子
                // for (int neighbor_cell : _neighbor_cell_Array[cell_index])
                for (int ii = 0; ii < ni; ++ii)
                {
                    int nebcout = 0;
                    for (int neb = 0; neb < 27; neb++)
                    {
                        int neighbor_cell = _neighbor_cell_Array[cell_index][neb];
                        // std::cout << "neighbor_cell: " << neighbor_cell << std::endl/;
                        std::array<int,3> shift = _neighbor_cell_Array_shifts[cell_index][neb];
                        
                        // for (int j : _cell_list.get_atoms_in_cell(neighbor_cell))        // Loop through the atoms in the neighboring cells
                        int nj = _count[neighbor_cell];                        // 查找单元内有几个原子
                        int j = _head[neighbor_cell];                         // 根据单元序号查找单元的第一个原子
                        for (int jj = 0; jj < nj; ++jj)
                        {
                            std::array<float,3> pos_i = xyz[i];
                            std::array<float,3> pos_j = xyz[j];
                            float r2 = _box.calc_distance2(pos_i, pos_j);
                            // Vec3<float> dr = pos_i - pos_j - shift * box_length;
                            // float r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
                            if ((r2 - rc2 < 1e-7) && (i != j))
                            {
                                // std::cout << "r2: " << r2 << "  " << j << std::endl;
                                // _neighlist_1d[_nnebs * i + nebcout] = j;
                                // nebcout++;
                                _neighborListArray[i].push_back(j);
                            }
                            j = _lscl[j];
                        }
                    }
                    i = _lscl[i];
                }
            }
            
            auto build_end = std::chrono::high_resolution_clock::now();
            auto build_duration = std::chrono::duration_cast<std::chrono::milliseconds>(build_end - build_start);
            //std::cout << "build cell list time: " << build_duration.count() << "ms" << std::endl;
        }
    }


    std::vector<std::vector<int>> CellSearchCPU::get_particle_neighbor_info()
    {
        // for (int i = 0; i < _natoms; i++)
        // {
        //     std::vector<int> neighbors;
        //     for (int j = 0; j < _nnebs; j++)
        //     {
        //         if (_neighlist_1d[_nnebs * i + j] != -1)
        //         {
        //             neighbors.push_back(_neighlist_1d[_nnebs * i + j]);
        //         }
        //     }
        //     _neighborListArray[i] = neighbors;
        // }
        return _neighborListArray;
    }


} // namespace dpnblist

