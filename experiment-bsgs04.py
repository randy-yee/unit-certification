read("ExperimentFunctions.py");
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ Starting precision:
{
sigstring = "0-4";



b_ranges =
[
    1,6,1;          \\#best 4 very small sized reg (~7)
    1,6,1;          \\#best 3
    1,6,1;          \\#best 3
    5,30,2;          \\#best 9 timeout 50
    2,25,2;          \\#best 6 timeout 50
    4,25,2;          \\#best 5 timeout 50
    6,25,2;          \\#best 12 timeout 50
    4,25,2;          \\#best 14 timeout 50
    6,20,2;          \\#best 14  timeout 50
    40, 80, 4;      \\#best 60, timeout 200
    40, 90, 5;      \\#best 45, timeout 200
    10, 60, 5;      \\#best 20, timeout 100
    200, 400, 50;     \\#best 300, timeout 500
    100, 200, 10;     \\#best 150, timeout 500
    80, 180, 20     \\#best 100, timeout 500
];
}
