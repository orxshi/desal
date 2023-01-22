#include <algorithm>
#include "patch.h"

namespace Liver
{
    Patch::Patch(PatchName name): name_(name)
    {
    }

    PatchName Patch::name() const
    {
        return name_;
    }

    void Patch::add_cell(Cell cell, NumberOfCells ncells)
    {
        if (cells_.capacity() < ncells)
        {
            cells_.reserve(ncells);
        }

        cells_.push_back(cell);
    }

    void Patches::add_cell(Cell cell, NumberOfCells ncells)
    {
        PatchName name = cell.patch_;

        auto patch = std::find_if(patches_.begin(), patches_.end(), [&](const Patch& f){return f.name() == name;});

        if (patch == patches_.end())
        {
            patches_.push_back(Patch(name));
            patches_.back().add_cell(cell, ncells);
        }
        else
        {
            patch->add_cell(cell, ncells);
        }
    }

    Patch* Patches::patch(PatchName name)
    {
        auto patch = std::find_if(patches_.begin(), patches_.end(), [&](const Patch& f){return f.name() == name;});

        if (patch == patches_.end())
        {
            return nullptr;
        }

        return &(*patch);
    }

    Cells& Patch::cells()
    {
        return cells_;
    }

    std::vector<Patch>::iterator Patches::begin()
    {
        return patches_.begin();
    }

    std::vector<Patch>::iterator Patches::end()
    {
        return patches_.end();
    }
}
