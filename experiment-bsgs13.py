read("ExperimentFunctions.py");
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ Starting precision:
{
sigstring = "1-3";



b_ranges =
[
    1,5, 1;                      \\# [10, 40] best 2
    1,7, 1;                      \\# [10, 40] best 5
    1,7, 1;                      \\# [10, 40] best 2
    1,5, 1;                      \\# [10, 40] best 2
    1,7, 1;                      \\# [10, 40] best 2 (these fields are too small)
    15,30, 2;                   \\# [10, 500] best 21
    15,25, 2;                   \\# [10, 500] best 19
    2,12, 2;                    \\# [10, 500] best 6
    2,12, 2;                    \\# [10, 500] best 10
    20,30, 2;                   \\# [10, 500] best 26
    150,  160, 5;                 \\#[1000, 500] best 160
    170, 195, 5;                  \\#[1000, 500] best 190
    110, 140, 5;                  \\#[1000, 500] best 130
    65,  85, 5;                   \\#[1000, 500] best 75
    120, 160, 5;                  \\#[1000, 500] best 125
    500, 1000, 100              \\
];

b_ranges_round2 =
[
    1,10, 1;                      \\# [10, 40] best 2
    1,7, 1;                      \\# [10, 40] best 5
    1,10, 1;                      \\# [10, 40] best 2
    1,7, 1;                      \\# [10, 40] best 2
    1,6, 1;                      \\# [10, 40] best 2 (these fields are too small)
    50,80, 5;                   \\# [10, 500] best 21
    50,75, 5;                   \\# [10, 500] best 55
    25,45, 5;                   \\# [10, 500] best 6
    40,60, 5;                   \\# [10, 500] best 10
    25,45, 5;                   \\# [10, 500] best 26
    450, 550, 25;                 \\#[1000, 500] best 160
    300, 400, 25;                  \\#[1000, 500] best 190
    450, 550, 25;                  \\#[1000, 500] best 130
    550, 650, 25;                  \\#[1000, 500] best 75
    375, 475, 25                   \\#[1000, 500] best 125
];
}
