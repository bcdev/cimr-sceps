package org.esa.cimr.sceps;

import java.io.File;

/**
 * Constants for cimr.sceps module
 *
 * @author olafd
 *
 */
public class ScepsConstants {

    public static final String SCEPS_CONFIG_ELEMENTS_TAG_NAME = "parameter";
    public static final String SCEPS_CODES_FOLDER_NAME = "SCEPScodes";
    public static final String SCEPS_CODES_OSFI_MATLAB_FOLDER_NAME =
             "OSFI" + File.separator + "Matlab";
    public static final String SCEPS_CODES_GENERAL_SUBMODULES_FOLDER_NAME =
            "General" + File.separator + "SubModules";
    public static final String SCEPS_DATA_FOLDER_NAME = "SCEPSdata";

    public static final String SCENE_GENERATION_MODULE_FOLDER_NAME = "SceGenMod";
    public static final String SCENE_GENERATION_MODULE_MODULES_SUBFOLDER_NAME = "Modules";
    public static final String SCENE_GENERATION_MODULE_SUBMODULES_SUBFOLDER_NAME = "SubModules";
    public static final String DEVALGO_L2_MODULE_FOLDER_NAME = "Devalgo_L2";
    public static final String DEVALGO_L2_MODULE_MODULES_SUBFOLDER_NAME = "Modules";

    public static final String SCEPS_SCD_ROOT_CONFIG_ITEM_NAME = "sceps_scd_root";
    public static final String SCEPS_MODULE_NAME_CONFIG_ITEM_NAME = "module_name";
    public static final String SCEPS_SCENE_TYPE_CONFIG_ITEM_NAME = "scene_type";
    public static final String SCEPS_SCENE_DATE_CONFIG_ITEM_NAME = "scene_date";

    public static final String SCEPS_INPUT_DATA_FOLDER_CONFIG_ITEM_NAME = "input_data_folder";
    public static final String SCEPS_INPUT_DATA_FILENAME_CONFIG_ITEM_NAME = "input_data_filename";
    public static final String SCEPS_OUTPUT_DATA_FOLDER_CONFIG_ITEM_NAME = "output_data_folder";
    public static final String SCEPS_OUTPUT_DATA_FILENAME_CONFIG_ITEM_NAME = "output_data_filename";

    public static final String SCEPS_DEVALGO_L2_L2GRID_CONFIG_ITEM_NAME = "l2_grid";

}
