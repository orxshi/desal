#ifndef PATCH_H
#define PATCH_H

#include "cell.h"

namespace Liver
{
    class Patch
    {
        public:

            Patch(PatchName);

            void add_cell(Cell, NumberOfCells);
            PatchName name() const;
            Cells& cells();

            PatchName name_;
            Cells cells_;

        private:

    };

    class Patches
    {
        public:

            void add_cell(Cell, NumberOfCells);
            Patch* patch(PatchName);
            std::vector<Patch>::iterator begin();
            std::vector<Patch>::iterator end();

            std::vector<Patch> patches_;

        private:

    };
}

#endif
