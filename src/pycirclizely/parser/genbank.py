from __future__ import annotations

import bz2
import gzip
import warnings
import zipfile
from collections import defaultdict
from io import StringIO, TextIOWrapper
from pathlib import Path
from typing import TYPE_CHECKING, Optional, Union

import numpy as np
from Bio import SeqIO, SeqUtils
from Bio.SeqFeature import Seq, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

if TYPE_CHECKING:
    from numpy.typing import NDArray


class Genbank:
    """Genbank Parser Class"""

    def __init__(
        self,
        gbk_source: str | Path | TextIOWrapper | list[SeqRecord],
        *,
        name: str | None = None,
        min_range: None = None,
        max_range: None = None,
    ):
        """Initialize Genbank.

        Args:
            gbk_source: Genbank file or source
                (`*.gz`, `*.bz2`, `*.zip` compressed file can be readable).
            name: Name (If None, `file name` or `record name` is set).
            min_range: No longer used. Left for backward compatibility.
            max_range: No longer used. Left for backward compatibility.
        """
        self._gbk_source = gbk_source
        if isinstance(gbk_source, (str, Path, StringIO, TextIOWrapper)):
            self._records = self._parse_gbk_source(gbk_source)
        else:
            self._records = gbk_source

        # Set genbank name
        if name is not None:
            self._name = name
        elif isinstance(self._gbk_source, (str, Path)):
            gbk_file = Path(self._gbk_source)
            if gbk_file.suffix in (".gz", ".bz2", ".zip"):
                self._name = gbk_file.with_suffix("").with_suffix("").name
            else:
                self._name = gbk_file.with_suffix("").name
        elif isinstance(self._gbk_source, (StringIO, TextIOWrapper)):
            self._name = self._records[0].name
        else:
            raise ValueError("Failed to get genbank name.")

        if min_range or max_range:
            warnings.warn("min_range & max_range is no longer used in Genbank parser.")

        if len(self.records) == 0:
            raise ValueError(f"Failed to parse {gbk_source} as Genbank file.")

    ############################################################
    # Property
    ############################################################

    @property
    def name(self) -> str:
        """Name"""
        return self._name

    @property
    def records(self) -> list[SeqRecord]:
        """Genbank records"""
        return self._records

    @property
    def genome_seq(self) -> str:
        """Genome sequence (only first record)"""
        return str(self.records[0].seq)

    @property
    def genome_length(self) -> int:
        """Genome length (only first record)"""
        return len(self.genome_seq)

    @property
    def range_size(self) -> int:
        """Same as `self.genome_length` (Left for backward compatibility)"""
        return self.genome_length

    @property
    def full_genome_seq(self) -> str:
        """Full genome sequence (concatenate all records)"""
        return "".join(str(r.seq) for r in self.records)

    @property
    def full_genome_length(self) -> int:
        """Full genome length (concatenate all records)"""
        return len(self.full_genome_seq)

    ############################################################
    # Public Method
    ############################################################

    def calc_genome_gc_content(self, seq: str | None = None) -> float:
        """
        Args:
            seq: Sequence for GC content calculation (Default: `self.genome_seq`).
        """
        seq = self.genome_seq if seq is None else seq
        gc_content = SeqUtils.gc_fraction(seq) * 100
        return gc_content

    def calc_gc_skew(
        self,
        window_size: Optional[int] = None,
        step_size: Optional[int] = None,
        *,
        seq: Optional[str] = None,
    ) -> tuple[NDArray[np.int64], NDArray[np.float64]]:
        """
        Args:
            window_size: Window size (Default: `genome_size / 500`).
            step_size: Step size (Default: `genome_size / 1000`).
            seq: Sequence for GCskew calculation (Default: `self.genome_seq`).
        """
        seq = self.genome_seq if seq is None else seq
        if window_size is None:
            window_size = int(len(seq) / 500)
        if step_size is None:
            step_size = int(len(seq) / 1000)
        if window_size == 0 or step_size == 0:
            window_size, step_size = len(seq), int(len(seq) / 2)
        pos_list = np.arange(0, len(seq), step_size, dtype=np.int64)
        if len(pos_list) == 0 or pos_list[-1] != len(seq):
            pos_list = np.append(pos_list, len(seq))

        # Initialize GC skew array and compute
        gc_skew_list = np.empty(len(pos_list), dtype=np.float64)
        for i, pos in enumerate(pos_list):
            start = max(0, pos - window_size // 2)
            end = min(len(seq), pos + window_size // 2)
            subseq = seq[start:end]

            g = subseq.count("G") + subseq.count("g")
            c = subseq.count("C") + subseq.count("c")

            try:
                skew = (g - c) / float(g + c)
            except ZeroDivisionError:
                skew = 0.0

            gc_skew_list[i] = skew

        return pos_list, gc_skew_list

    def calc_gc_content(
        self,
        window_size: Optional[int] = None,
        step_size: Optional[int] = None,
        *,
        seq: Optional[str] = None,
    ) -> tuple[NDArray[np.int64], NDArray[np.float64]]:
        """
        Args:
            window_size: Window size (Default: `genome_size / 500`).
            step_size: Step size (Default: `genome_size / 1000`).
            seq: Sequence for GC content calculation (Default: `self.genome_seq`).
        """
        seq = self.genome_seq if seq is None else seq
        assert seq is not None

        # Handle default values
        if window_size is None:
            window_size = int(len(seq) / 500)
        if step_size is None:
            step_size = int(len(seq) / 1000)
        if window_size == 0 or step_size == 0:
            window_size, step_size = len(seq), int(len(seq) / 2)

        positions = list(range(0, len(seq), step_size))
        if not positions or positions[-1] != len(seq):
            positions.append(len(seq))
        pos_list = np.array(positions, dtype=np.int64)

        # Initialize and fill GC content array
        gc_content_list = np.empty_like(pos_list, dtype=np.float64)
        for i, pos in enumerate(pos_list):
            window_start_pos = pos - int(window_size / 2)
            window_end_pos = pos + int(window_size / 2)
            window_start_pos = max(0, window_start_pos)
            window_end_pos = min(len(seq), window_end_pos)

            subseq = seq[window_start_pos:window_end_pos]
            gc_content_list[i] = SeqUtils.gc_fraction(subseq) * 100

        return pos_list, gc_content_list

    def get_seqid2seq(self) -> dict[str, str]:
        """Get seqid & complete/contig/scaffold genome sequence dict."""
        return {str(rec.id): str(rec.seq) for rec in self.records}

    def get_seqid2size(self) -> dict[str, int]:
        """Get seqid & complete/contig/scaffold genome size dict."""
        return {seqid: len(seq) for seqid, seq in self.get_seqid2seq().items()}

    def get_seqid2features(
        self,
        feature_type: Union[str, list[str], None] = "CDS",
        target_strand: Optional[int] = None,
    ) -> dict[str, list[SeqFeature]]:
        """Get seqid & features in target seqid genome dict

        Args:
            feature_type: Feature type (`CDS`, `gene`, `mRNA`, etc...)
                If None, extract regardless of feature type.
            target_strand: Extract target strand.
                If None, extract regardless of strand.
        """
        if isinstance(feature_type, str):
            feature_type = [feature_type]

        seqid2features: dict[str, list[SeqFeature]] = defaultdict(list)
        for rec in self.records:
            for feature in rec.features:
                # Ignore feature if parsing of location fails
                if feature.location is None:
                    continue

                # Filter feature by type & strand
                strand = feature.location.strand
                if feature_type is not None and feature.type not in feature_type:
                    continue
                if target_strand is not None and strand != target_strand:
                    continue
                # Exclude feature which straddle genome start position
                if self._is_straddle_feature(feature):
                    continue

                start = int(feature.location.start)
                end = int(feature.location.end)
                seqid2features[str(rec.id)].append(
                    SeqFeature(
                        location=SimpleLocation(start, end, strand),
                        type=feature.type,
                        qualifiers=feature.qualifiers,
                    )
                )
        return dict(seqid2features)  # Convert defaultdict to dict

    def extract_features(
        self,
        feature_type: str | list[str] | None = "CDS",
        *,
        target_strand: int | None = None,
        target_range: tuple[int, int] | None = None,
    ) -> list[SeqFeature]:
        """Extract features (only first record)

        Args:
            feature_type: Feature type to extract.
            target_strand: Strand to extract. If None, extracts regardless of strand.
            target_range: Range to extract. If None, extracts regardless of range.
        """
        seqid2features = self.get_seqid2features(feature_type, target_strand)
        first_record_features = list(seqid2features.values())[0]
        if target_range:
            target_features = []
            for feature in first_record_features:
                start = int(feature.location.start)  # type: ignore
                end = int(feature.location.end)  # type: ignore
                if min(target_range) <= start <= end <= max(target_range):
                    target_features.append(feature)
            return target_features
        else:
            return first_record_features

    def write_cds_fasta(self, outfile: str | Path) -> None:
        """
        Args:
            outfile: Output CDS fasta file.
        """
        cds_records: list[SeqRecord] = []
        counter = 0
        seqid2cds_features = self.get_seqid2features(feature_type="CDS")
        for seqid, cds_features in seqid2cds_features.items():
            for cds_feature in cds_features:
                # Ignore no translation feature
                translation = cds_feature.qualifiers.get("translation", [None])[0]
                if translation is None:
                    continue
                # Get feature location
                start = int(cds_feature.location.start)  # type: ignore
                end = int(cds_feature.location.end)  # type: ignore
                strand = -1 if cds_feature.location.strand == -1 else 1
                # Set feature id
                location_id = f"|{seqid}|{start}_{end}_{strand}|"
                protein_id = cds_feature.qualifiers.get("protein_id", [None])[0]
                if protein_id is None:
                    feature_id = f"GENE{counter:06d}{location_id}"
                else:
                    feature_id = f"GENE{counter:06d}_{protein_id}{location_id}"
                counter += 1
                # Add SeqRecord of CDS feature
                seq = Seq(translation)
                product = cds_feature.qualifiers.get("product", [""])[0]
                seq_record = SeqRecord(seq, feature_id, description=product)
                cds_records.append(seq_record)
        # Write CDS file
        SeqIO.write(cds_records, outfile, "fasta-2line")

    def write_genome_fasta(self, outfile: str | Path) -> None:
        """
        Args:
            outfile: Output genome fasta file.
        """
        with open(outfile, "w", encoding="utf-8") as f:
            for seqid, seq in self.get_seqid2seq().items():
                f.write(f">{seqid}\n{seq}\n")

    ############################################################
    # Private Method
    ############################################################

    def _parse_gbk_source(
        self, gbk_source: str | Path | TextIOWrapper
    ) -> list[SeqRecord]:
        """Parse genbank source."""
        if isinstance(gbk_source, (str, Path)):
            if Path(gbk_source).suffix == ".gz":
                with gzip.open(gbk_source, mode="rt", encoding="utf-8") as f:
                    return list(SeqIO.parse(f, "genbank"))
            elif Path(gbk_source).suffix == ".bz2":
                with bz2.open(gbk_source, mode="rt", encoding="utf-8") as f:
                    return list(SeqIO.parse(f, "genbank"))
            elif Path(gbk_source).suffix == ".zip":
                with zipfile.ZipFile(gbk_source) as zip:
                    with zip.open(zip.namelist()[0]) as f:
                        io = TextIOWrapper(f, encoding="utf-8")
                        return list(SeqIO.parse(io, "genbank"))
            else:
                with open(gbk_source, encoding="utf-8") as f:
                    return list(SeqIO.parse(f, "genbank"))
        # Parse TextIOWrapper
        return list(SeqIO.parse(gbk_source, "genbank"))

    def _is_straddle_feature(self, feature: SeqFeature) -> bool:
        """Check target feature straddle genome start position or not."""
        strand = feature.location.strand
        if strand == -1:
            start = int(feature.location.parts[-1].start)  # type: ignore
            end = int(feature.location.parts[0].end)  # type: ignore
        else:
            start = int(feature.location.parts[0].start)  # type: ignore
            end = int(feature.location.parts[-1].end)  # type: ignore
        return True if start > end else False

    def __str__(self):
        text = f"{self.name}: {len(self.records)} records\n"
        for num, (seqid, size) in enumerate(self.get_seqid2size().items(), 1):
            text += f"{num:02d}. {seqid} ({size:,} bp)\n"
        return text
