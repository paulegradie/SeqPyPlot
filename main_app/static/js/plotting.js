$('textarea').keypress(function(event) {
    if (event.which == 13) {
        event.stopPropagation();
    }
});

$('#plotGeneList').keypress(function(event) {
    event.stopPropagation();
})

$(function(){
    $("#select").change(function(){
        var text = "";
        $( "select option:selected" ).each(function() {
            text += $(this).text() += " ";
        })
        $( "#plottableInputText" ).text( text )
    }).change();
})

function plotting_form_validation (curForm) {

    function validate_extensions (curForm) {
        var valid_extensions = ['spp', 'txt']
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
    function validate_gene_list (curForm) {
        var gene_list_inputs = document.getElementsByName("gene_list")
        console.log(gene_list_inputs)
        content = gene_list_inputs[0].value
        var genelist_not_empty = false;
        if (content.length > 0) {
            genelist_not_empty = true;
        }
        return genelist_not_empty
    }

    function validate_file_present (curForm, element_name) {
        var form_input = document.getElementsByName(element_name);

        var file_is_present = false;
        var filename = form_input[0].value
        if (filename.length > 0) {
            file_is_present = true
        }
        return file_is_present
    }

    // calls
    var data_file_present = validate_file_present(curForm, "data_file")
    var analysis_config_present = validate_file_present(curForm, "analysis_config")
    var DE_genelist_present = validate_file_present(curForm, "de_gene_list")
    var gene_list_present = validate_gene_list(curForm)
    var extensions_valid = validate_extensions(curForm)

    var form_validated = true;
    if(!data_file_present) {
        alert("Plottable data file is missing.")
        form_validated = false;
    }
    if (!analysis_config_present) {
        alert("Config file is missing.")
        form_validated = false;
    }
    if (!gene_list_present) {
        alert("Gene list is emtpy.")
        form_validated = false;
    }
    if(!extensions_valid) {
        alert('Only files with extension .spp allowed')
        form_valiated = false;
    }
    return true; //form_validated;
}
