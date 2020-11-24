#!/usr/bin/env python

"""
.. See the NOTICE file distributed with this work for additional information
   regarding copyright ownership.

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
"""
import argparse
import sys

from basic_modules.workflow import Workflow
from utils import logger
from apps.jsonapp import JSONApp
from tool.VRE_Tool import RUNNER


class process_WF_RUNNER(Workflow):
    """
    Functions for demonstrating the pipeline set up.
    """

    configuration = {}
    output = {}

    def __init__(self, configuration=None):
        """
        Initialise the tool with its configuration.

        :param configuration: a dictionary containing parameters that define how the operation should be carried out,
        which are specific to each Dorothea.
        :type configuration: dict
        """
        logger.debug("Processing Dorothea")
        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    def run(self, input_files, metadata, output_files, output_metadata):
        """
        Main run function for processing a test file.

        :param input_files: Dictionary of file locations.
        :type input_files: dict
        :param metadata: Required meta data.
        :type metadata: dict
        :param output_files: Locations of the output files to be returned by the pipeline.
        :type output_files: dict
        :param output_metadata:
        :type output_metadata: list
        :return: Locations for the output txt (output_files), Matching metadata for each of the files (output_metadata).
        :rtype: dict, dict
        """
        try:
            logger.debug("Initialise the Dorothea")
            tt_handle = RUNNER(self.configuration)
            tt_files, tt_meta = tt_handle.run(input_files, metadata, output_files, output_metadata)
            return tt_files, tt_meta

        except Exception as error:
            errstr = "Test Dorothea wasn't processed successfully. ERROR: {}".format(error)
            logger.error(errstr)
            raise Exception(errstr)


def main_json(config_path, in_metadata_path, out_metadata_path):
    """
    Main function.

    This function launches the app using configuration written in two json files: defineIO.json and input_metadata.json.

    :param config_path: path to a valid JSON file containing information on how the tool should be executed.
    :type config_path: str
    :param in_metadata_path: path to a valid JSON file containing information on tool inputs.
    :type in_metadata_path: str
    :param out_metadata_path: path to write the JSON file containing information on tool outputs.
    :type out_metadata_path: str
    :return: If result is True, execution finished successfully. False, otherwise.
    :rtype: bool
    """
    try:
        logger.info("1. Instantiate and launch the App")
        app = JSONApp()

        result = app.launch(process_WF_RUNNER, config_path, in_metadata_path, out_metadata_path)  # launch the app
        logger.info("2. App successfully launched; see " + out_metadata_path)
        return result

    except Exception as error:
        errstr = "App wasn't successfully launched. ERROR: {}".format(error)
        logger.error(errstr)
        raise Exception(errstr)


if __name__ == "__main__":

    # Set up the command line parameters
    PARSER = argparse.ArgumentParser(description="VRE Dorothea runner")
    PARSER.add_argument("--config", help="Configuration file", required=True)
    PARSER.add_argument("--in_metadata", help="Location of input metadata file", required=True)
    PARSER.add_argument("--out_metadata", help="Location of output metadata file", required=True)
    PARSER.add_argument("--log_file", help="Location of the log file", required=False)

    # Get the matching parameters from the command line
    ARGS = PARSER.parse_args()

    CONFIG = ARGS.config
    IN_METADATA = ARGS.in_metadata
    OUT_METADATA = ARGS.out_metadata

    if ARGS.log_file:
        sys.stderr = sys.stdout = open(ARGS.log_file, "w")

    RESULTS = main_json(CONFIG, IN_METADATA, OUT_METADATA)
