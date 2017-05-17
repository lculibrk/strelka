//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2017 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

///
/// \author Sangtae Kim
///


#include "CandidateSnvBuffer.hh"

void CandidateSnvBuffer::addCandidateSnv(
    const unsigned sampleId,
    const pos_t pos,
    const char baseChar,
    const HaplotypeId haplotypeId,
    const float altHaplotypeCountRatio)
{
    assert (haplotypeId == 1 || haplotypeId == 2);
    HaplotypeIdAndCountRatio& haplotypeIdForBase = _candidateSnvBuffer[sampleId].getRef(pos);
    switch (baseChar)
    {
    case 'A':
        haplotypeIdForBase.a += haplotypeId;
        break;
    case 'C':
        haplotypeIdForBase.c += haplotypeId;
        break;
    case 'G':
        haplotypeIdForBase.g += haplotypeId;
        break;
    case 'T':
        haplotypeIdForBase.t += haplotypeId;
        break;
    default:
        assert(false);
    }
    haplotypeIdForBase.altHaplotypeCountRatio += altHaplotypeCountRatio;
}

bool CandidateSnvBuffer::isCandidateSnv(const unsigned sampleId, const pos_t pos, const char baseChar) const
{
    const auto baseIndex(static_cast<BASE_ID::index_t>(base_to_id(baseChar)));
    if (baseIndex == BASE_ID::ANY) return false;
    return (getHaplotypeId(sampleId, pos, baseIndex) != 0);
}

HaplotypeId CandidateSnvBuffer::getHaplotypeId(const unsigned sampleId, const pos_t pos, const BASE_ID::index_t baseIndex) const
{
    if (not _candidateSnvBuffer[sampleId].isKeyPresent(pos))
    {
        return static_cast<HaplotypeId>(0); // haplotype ID 0 means there's no candidate SNV at the position
    }

    const auto haplotypeIdForBase(_candidateSnvBuffer[sampleId].getConstRef(pos));
    switch (baseIndex)
    {
    case BASE_ID::A:
        return haplotypeIdForBase.a;
    case BASE_ID::C:
        return haplotypeIdForBase.c;
    case BASE_ID::G:
        return haplotypeIdForBase.g;
    case BASE_ID::T:
        return haplotypeIdForBase.t;
    default:
        assert(false);
    }
    return static_cast<HaplotypeId>(0);
}

float CandidateSnvBuffer::getAltHaplotypeCountRatio(const unsigned sampleId, const pos_t pos) const
{
    if (not _candidateSnvBuffer[sampleId].isKeyPresent(pos))
    {
        return 0.0f;
    }
    return _candidateSnvBuffer[sampleId].getConstRef(pos).altHaplotypeCountRatio;
}

bool CandidateSnvBuffer::empty() const
{
    return _candidateSnvBuffer.empty();
}

void CandidateSnvBuffer::clearUpToPos(const pos_t pos)
{
    for (unsigned sampleId(0); sampleId<_candidateSnvBuffer.size(); ++sampleId)
    {
        _candidateSnvBuffer[sampleId].eraseTo(pos);
    }
}