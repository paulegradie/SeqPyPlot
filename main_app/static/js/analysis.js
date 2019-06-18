    function get(key, map) {
        return map[key]
    }

    $("#singlesample-search").hide();
    let defaultConditions = 2
    window.current_form = "#doublesample-search"

    let single_file_rows = ["singlecontrol_1"];
    let single_file_names = ["singletimename_1"];

    let double_file_rows_control = ["doublecontrol_1"];
    let double_file_rows_treated = ["doubletrated_1"];
    let double_file_names = ["doubletimename_1"];

    let single_params = ["param1_1", "param1_2", "param1_3", "param1_4", "param1_5", "param1_6", "param1_7", "param1_8"];
    let double_params = ["param2_1", "param2_2", "param2_3", "param2_4", "param2_5", "param2_6", "param2_7", "param2_8"];


    $("select#numCondition").change(function () {
        let map = {'one': 1, 'two': 2}
        function get (map, key) {
            return map[key];
        }

        var numconditions_str = $(this).children("option:selected").val()
        // var numconditions_int = map[numconditions]
        var numconditions_int = get(map, numconditions_str)

        if (numconditions_str == "one") {
            console.log('OK')
            $("#doublesample-search").hide().attr("formnovalidate");
            $("#singlesample-search").toggle();
            window.current_form = "#singlesample-search"
            console.log(window.current_form)

        } else if (numconditions_str == "two") {
            $("#singlesample-search").hide().attr("formnovalidate");
            $("#doublesample-search").toggle();
            window.current_form = "#doublesample-search"
            console.log(window.current_form)

        } else {
            alert("MUST CHOOSE")
        }
    })

    // var handles the numebr of sample rows
    let curSample = 2

    // add a row to the sample table
    $("#addSample").click(function () {

        var x1 = $(this).offset().top;

        $(function () {
            if (curSample <= 12) {

                var rowNumber = curSample.toString()

                var newRowSingle='<tr><td><span class="sampleItemSpan">'+rowNumber+'.</span></td><td><input name="time_name1_'+rowNumber+'" id="singletimename_'+rowNumber+'" class="sampleItemSpan" type="text" placeholder="T'+rowNumber+'"></td><td><input id="singlecontrol_'+rowNumber+'" name="control_'+rowNumber+'" class="sampleItemSpan" type="file"></td></tr>'
                var newRowDouble = '<tr><td><span class="">'+rowNumber+'.</span></td><td><input name="time_name2_'+rowNumber+'" id="doubletimename_'+rowNumber+'" class="sampleItemSpan" type="text" placeholder="T'+rowNumber+'"></td><td><input name="doublecontrol_'+rowNumber+'" class="sampleItemSpan" type="file"></td><td><input name="doubletreated_'+rowNumber+'" type="file"></td></tr>'

                $("#OneSampleTable").append(newRowSingle)
                $("#TwoSampleTable").append(newRowDouble)

                single_file_rows.push("singlecontrol_"+rowNumber);
                single_file_names.push("singletimename_"+rowNumber);

                double_file_rows_control.push("doublecontrol_"+rowNumber);
                double_file_rows_treated.push("doubletrated_"+rowNumber);
                double_file_names.push("doubletimename_"+rowNumber);
                console.log(double_file_names)
                curSample += 1

            } else {
                alert("Plots look bad with too many samples...")

            }
        })

        var x2 = $(this).offset().top;
        var dx = x2 - x1;
        $(document).scrollTop($(document).scrollTop() + dx);
    })

    $("#deleteSample").click(function () {

        var x1 = $(this).offset().top;
        $(function () {
            if (curSample > 2) {
                $("#OneSampleTable tr:last").remove()
                $("#TwoSampleTable tr:last").remove()
                curSample -= 1

                single_file_rows.pop()
                single_file_names.pop()

                double_file_rows_control.pop()
                double_file_rows_treated.pop()
                double_file_names.pop()
                console.log(double_file_names)


            } else {
                alert("Need at least 1 time point")
            }
            return false
        })
        var x2 = $(this).offset().top;
        var dx = x2 - x1;
        $(document).scrollTop($(document).scrollTop() + dx);
    })

    // main validation function
    function analysis_form_validation (curForm) {
        let valid_extensions = ['counts'];

        function validate_texts (text_ids) {
            var texts_present = true;

            for (var i=0; i<text_ids.length; i++) {
                cur_id = text_ids[i];
                el = document.getElementById(cur_id)
                if (el.value.length == 0) {
                    texts_present = false;
                }
            }
            return texts_present;
        }

        function validate_all_files_there (curForm) {

            form_inputs = curForm.getElementsByTagName("Input");
            var num_files_found = 0
            let file_count_is_good = false;

            if (window.current_form == '#doublesample-search') {
                var cur_cond_count = 2;
            } else if (window.current_form == '#singlesample-search') {
                var cur_cond_count = 1
            }

            for (var i = 0; i < form_inputs.length; i++) {
                var form_input = form_inputs[i];
                if (form_input.type == "file") {
                    var filename = form_input.value
                    if (filename.length > 0) {
                        num_files_found += 1
                    }
                }
            }
            if (num_files_found % cur_cond_count == 0) {
                    file_count_is_good = true;
                }
            if (num_files_found == 0) {
                file_count_is_good = false;
            }

            return file_count_is_good
        }

        function validate_extensions (curForm) {

            var form_inputs = curForm.getElementsByTagName("Input");
            var extensions_all_good = true;

            // get all extensions
            var extensions = [];
            for (var i = 0; i < form_inputs.length; i++) {
                var form_input = form_inputs[i]
                if (form_input.type == "file") {
                    filename = form_input.value
                    var filename_parts = filename.split('.')
                    var file_extension = filename_parts[filename_parts.length-1]
                    extensions.push(file_extension)
                }
            }

            // check each extension is in the allowed extensions
            for (var ext_idx = 0; ext_idx < extensions.length; ext_idx++) {

                var file_extension = extensions[ext_idx];
                if (valid_extensions.indexOf(file_extension) == -1) {
                    extensions_all_good = false
                }

            }
            return extensions_all_good
        }

        function validate_form_params (param_ids) { // param_ids is a list of param ids passed given the window.curren_form

            var all_params_present = true;

            for (var i = 0; i < param_ids.length; i++) {
                var param_id = param_ids[i];
                var param = document.getElementById(param_id)
                console.log("PARAM: " + param + ": " + param_id)

                if (param.value.length == 0) {
                    console.log("NO PARAM: "+param.value.length)
                    console.log(param.value)
                    all_params_present = false;
                }
            }
            return all_params_present;
        }

        // form_inputs = curForm.getElementsByTagName("Input");

        var all_files_present = validate_all_files_there(curForm);
        var extensions_valid = validate_extensions(curForm);

        if (window.current_form == '#doublesample-search') {
            var all_numbers_present = validate_form_params(double_params);
            var all_texts_present = validate_texts(double_file_names)

        } else if (window.current_form == '#singlesample-search') {
            var all_numbers_present = validate_form_params(single_params);
            var all_texts_present = validate_texts(single_file_names)
        } else {
            alert("Only one or two conditions supported.")
        }

        var form_validated = true;
        if (!all_files_present) {
            alert("Files are missing. You may pass empty files when imputing a time point, but all file slots must be filled.")
            form_validated = false;
        }

        if (!extensions_valid) {
            alert("Files must be of type .counts. Row format is two column: text-tab-number-newline")
            form_validated = false;
        }

        if (!all_numbers_present) {
            alert("Filter params are missing.")
            form_validated = false;
        }

        if (!all_texts_present) {
            alert("The experiment and files must all be named.")
            form_validated = false;
        }

        return form_validated;
    }
