read("ExperimentFunctions.py");
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ Starting precision:
{
sigstring = "1-3";



b_ranges =
[
    1,15,1;          \\#best 1 timeout 10 (regulator ~2)
    1,10,1;          \\#best 2 timeout 10
    1,6,1;          \\#best 2  max 20
    2,12,1;          \\# best 8 max 50
    1,12,1;          \\# best 3 max 50 ceil(i/3)== 2
    1,10,1;          \\# best 2 max 50
    25,55,3;          \\#best 30 max 400 (also 50)
    5,25,2;          \\#best 17 max 400
    9,25,2;          \\#best 13 max 400
    75, 95, 2;      \\#best 95 max 1000
    30, 65, 5;      \\#best 40 max 1000
    60, 100, 5;      \\#best 80 max 1000
    100, 140, 5;     \\#best 120 max 1500
    160, 230, 10;     \\#best 170 sub 300 max 1500
    130, 160, 10     \\#best 120 max 1500
];
}
