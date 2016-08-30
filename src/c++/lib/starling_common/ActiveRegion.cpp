// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2016 Illumina, Inc.
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

#include <boost/algorithm/string.hpp>
#include "ActiveRegion.hh"

void ActiveRegion::insertHaplotypeBase(align_id_t alignId, pos_t pos, const std::string& base)
{
    if (!_alignIdToHaplotype.count(alignId))
    {
        // first occurrence of this alignment
        _alignIdToHaplotype[alignId] = std::string();
        for (int i=_start; i<pos; ++i)
            _alignIdToHaplotype[alignId] += missingPrefix;
    }
    _alignIdToHaplotype[alignId] += base;
    if (pos == _end)
        _alignIdReachingEnd.insert(alignId);
}

// decompose haplotypes into primitive alleles
void ActiveRegion::processHaplotypes(IndelBuffer& indelBuffer, RangeSet& polySites) const
{
    std::map<std::string, std::vector<align_id_t>> haplotypeToAlignIdSet;
    std::vector<std::pair<std::string, align_id_t>> softClippedReads;
    for (const auto& entry : _alignIdToHaplotype)
    {
        align_id_t alignId = entry.first;
        const std::string& haplotype(entry.second);

        // ignore if the read does not cover the start of the active region
        if (haplotype.empty() || haplotype[0] == missingPrefix) continue;

        // ignore if the read does not reach the end of the active region
        if (_alignIdReachingEnd.find(alignId) == _alignIdReachingEnd.end()) continue;

        if (_alignIdSoftClipped.find(alignId) != _alignIdSoftClipped.end())
        {
            softClippedReads.push_back(std::pair<std::string, align_id_t>(haplotype, alignId));
            continue;
        }

        if (!haplotypeToAlignIdSet.count(haplotype))
            haplotypeToAlignIdSet[haplotype] = std::vector<align_id_t>();
        haplotypeToAlignIdSet[haplotype].push_back(alignId);
    }

    // match soft-clipped reads to the haplotypes
    for (auto& entry : haplotypeToAlignIdSet)
    {
        const std::string& haplotype(entry.first);
        auto& alignIdList(entry.second);

        for (const auto& softClipEntry : softClippedReads)
        {
            const std::string& softClippedRead(softClipEntry.first);
            if (boost::starts_with(softClippedRead, haplotype)
                or boost::ends_with(softClippedRead, haplotype))
            {
                align_id_t alignId(softClipEntry.second);
                alignIdList.push_back(alignId);
            }
        }
    }

    // determine threshold to select 3 haplotypes with the largest counts
    unsigned largestCount = MinHaplotypeCount;
    unsigned secondLargestCount = MinHaplotypeCount;
    unsigned thirdLargestCount = MinHaplotypeCount;
    unsigned totalCount = 0;
    for (const auto& entry : haplotypeToAlignIdSet)
    {
        auto count = entry.second.size();

        totalCount += count;
        if (count > thirdLargestCount)
        {
            if (count > secondLargestCount)
            {
                if (count > largestCount)
                {
                    thirdLargestCount = secondLargestCount;
                    secondLargestCount = largestCount;
                    largestCount = (unsigned)count;
                }
                else
                {
                    thirdLargestCount = secondLargestCount;
                    secondLargestCount = (unsigned)count;
                }
            }
            else
                thirdLargestCount = (unsigned)count;
        }
    }

//    std::cout << '>' << _start+1 << '\t' << _end+1 << '\t' << _refSeq << std::endl;
    for (const auto& entry : haplotypeToAlignIdSet)
    {
        const std::string& haplotype(entry.first);
        if (haplotype.empty() || haplotype[0] == missingPrefix) continue;

        const auto& alignIdList(entry.second);
        auto count = alignIdList.size();

//        std::cout << haplotype << '\t' << count << std::endl;
        if (count >= thirdLargestCount and haplotype != _refSeq)
        {
            convertToPrimitiveAlleles(haplotype, alignIdList, totalCount,
                                      indelBuffer, polySites);
        }
    }
}

void ActiveRegion::convertToPrimitiveAlleles(
    const std::string& haploptypeSeq,
    const std::vector<align_id_t>& alignIdList,
    const unsigned totalReadCount,
    IndelBuffer& indelBuffer,
    RangeSet& polySites) const
{
    AlignmentResult<int> result;
    _aligner.align(haploptypeSeq.begin(),haploptypeSeq.end(),_refSeq.begin(),_refSeq.end(),result);
    const ALIGNPATH::path_t& alignPath = result.align.apath;

    pos_t referencePos = _start;
    pos_t haplotypePosOffset = 0;
    if (result.align.beginPos > 0)
    {
        assert(false && "Unexpected alignment segment");
    }

    std::vector<pos_t> mismatchPositions;
    unsigned numVariants(0);
    for (unsigned pathIndex(0); pathIndex<alignPath.size(); ++pathIndex)
    {
        const ALIGNPATH::path_segment& pathSegment(alignPath[pathIndex]);
        unsigned segmentLength = pathSegment.length;

        std::unique_ptr<IndelKey> indelKeyPtr = nullptr;
        switch (pathSegment.type)
        {
        case ALIGNPATH::SEQ_MATCH:
            referencePos += segmentLength;
            haplotypePosOffset += segmentLength;
            break;
        case ALIGNPATH::SEQ_MISMATCH:
            for (unsigned i(0); i<segmentLength; ++i)
            {
                mismatchPositions.push_back(referencePos);
                ++referencePos;
                ++haplotypePosOffset;
            }
            ++numVariants;
            break;
        case ALIGNPATH::INSERT:
        {
            auto insertSeq(haploptypeSeq.substr(haplotypePosOffset, segmentLength));
            indelKeyPtr = std::unique_ptr<IndelKey>(new IndelKey(referencePos, INDEL::INDEL, 0, insertSeq.c_str()));

            haplotypePosOffset += segmentLength;
            ++numVariants;
            break;
        }
        case ALIGNPATH::DELETE:
        {
            indelKeyPtr = std::unique_ptr<IndelKey>(new IndelKey(referencePos, INDEL::INDEL, segmentLength));
            referencePos += segmentLength;
            ++numVariants;
            break;
        }
        default:
            assert(false && "Unexpected alignment segment");
        }

        if (indelKeyPtr != nullptr)
        {
            auto* indelDataPtr = indelBuffer.getIndelDataPtr(*indelKeyPtr);
            if (indelDataPtr == nullptr)
            {
                // novel indel
                for (auto alignId : alignIdList)
                {
                    IndelObservationData indelObservationData;
                    auto& alignInfo(_alignIdToAlignInfo[alignId % MaxDepth]);
                    indelObservationData.iat = alignInfo.indelAlignType;
                    indelObservationData.id = alignId;
                    indelBuffer.addIndelObservation(alignInfo.sampleId, {*indelKeyPtr, indelObservationData});
                }
                indelDataPtr = indelBuffer.getIndelDataPtr(*indelKeyPtr);
            }
            assert(indelDataPtr != nullptr && "Missing indelData");

            // determine whether this indel is candidate or private
            indelDataPtr->isConfirmedInActiveRegion = true;
        }
    }

    unsigned numReads = alignIdList.size();
    if (numVariants == 1 and (numReads < totalReadCount*HaplotypeFrequencyThreshold))
        return;
    for (pos_t mismatchPos : mismatchPositions)
        polySites.getRef(mismatchPos) = 1;
}



