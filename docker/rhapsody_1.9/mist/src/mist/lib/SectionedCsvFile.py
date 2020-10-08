from typing import Sequence, Union, Dict
import csv


class SectionedCsvFile(object):

    def __init__(self, file_handle):
        self.file_handle = file_handle
        self.data_model = []

    def add_section(self, section_title: str, data: Sequence[dict]):
        """

        Args:
            section_title: the title of the csv subsection
            columns: the order the columns
            data: rows of the output

        Usage:
            >>> section_csv = SectionedCsvFile(open("foo.csv"))
            >>> section_csv.add_section("Read Metrics", [{"foo": 1, "bar": "bazz"}])

        """
        self.data_model.append((section_title, data))

    def write(self):
        for section_title, data in self.data_model:
            self.file_handle.write("#{}#\n".format(section_title))
            section_writer = csv.DictWriter(self.file_handle, fieldnames=data[0].keys())
            section_writer.writeheader()
            for row in data:
                section_writer.writerow(row)
            self.file_handle.write("\n")