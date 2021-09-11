
    $(document).ready(function () {

        //attach the click event to all <img> tag having class 'popos_img'
        //here no need of any for loop


        // Show / hide different parts of the form based on what type of file we're working with
        $('#type').change(function () {
            if ($(this).val() === 'species') {
                $("input[type='checkbox']").prop('disabled', false).prop('checked', true);
                $('.genome-search').css('color', 'black');
                $('#add_sequence').prop('disabled', true);
                $('#add_sequence').prop('checked', false);
                $("label[for='add_sequence']").css('color', '#ccc');

            }
            else if ($(this).val() === 'profile') {
                $("input[type='checkbox']").prop('disabled', true).prop('checked', false);
                $('.genome-search').css('color', '#ccc');
                $("label[for='add_sequence']").css('color', '#ccc');


            }

            else {
                $("input[type='checkbox']").prop('disabled', false).prop('checked', true);
                $('.genome-search').css('color', 'black');
                $("label[for='add_sequence']").css('color', 'black');



            }
        });





    });






