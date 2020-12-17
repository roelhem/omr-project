%% Test Scalar input
assert(strcmp(lib.utils.fixDateAutocorrect('0-Sep'), '0-9'), '"Sept" is converted incorrectly');
assert(strcmp(lib.utils.fixDateAutocorrect('Oct-19'), '10-19'), '"Oct" is converted incorrectly');
assert(strcmp(lib.utils.fixDateAutocorrect('20-29'), '20-29'), "Strings that were already formatted correctly don't return itself as it's result.");

%% Test array input.
assert(all(strcmp(lib.utils.fixDateAutocorrect(...
    [{'0-Sep'},{'Oct-19'},{'20-29'}]...
),[{'0-9'},{'10-19'},{'20-29'}]...
)), "Horizontal array's aren't converted correctly.")
assert(all(strcmp(lib.utils.fixDateAutocorrect([
    {'0-Sep'}
    {'Oct-19'}
    {'20-29'}
]),[
    {'0-9'}
    {'10-19'}
    {'20-29'}
])), "Vertical array's aren't converted correctly.")