
read("ExperimentFunctions.py");
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{

sigstring = "2-2";

b_ranges =
[
    2, 6, 1;                       \\# [1,  40] best 3, very bad 15+
    2, 6, 1;                       \\# [10, 40] best 5, very bad 15+
    1, 6, 1;                       \\# [10, 40] best 3
    1, 5, 1;                       \\# [10, 40] best 3,
    1, 6, 1;                       \\# [10, 40] best 2,
    25, 45, 2;                     \\# [10, 500] best 35
    44, 56, 2;                     \\# [10, 500] best 50
    45, 66, 2;                     \\# [10, 500] best 61
    28, 40, 2;                     \\# [10, 500] best 32
    26, 45, 2;                     \\# [10, 500] best 40
    500, 560, 10;                 \\#[1000, 2000] best 540
    230, 260, 10;                 \\#[1000, 2000] best 240
    580, 640, 10;                 \\#[1000, 2000] best 620
    730, 770, 10;                 \\#[1000, 2000] best 750
    220, 280, 10;                 \\#[1000, 2000] best 250
    1000, 2000, 200
];

b_ranges_round2 =
[
    5, 15, 2;                       \\# [1,  40] best 3, very bad 15+
    5, 15, 2;                       \\# [10, 40] best 5, very bad 15+
    5, 15, 2;                       \\# [10, 40] best 3
    3, 10, 2;                       \\# [10, 40] best 3,
    3, 9, 1;                        \\# [10, 40] best 2,
    110, 150, 5;                     \\# [10, 500] best 35
    80, 110, 5;                      \\# [10, 500] best 50
    120, 150, 5;                     \\# [10, 500] best 61
    120, 200, 10;                    \\# [10, 500] best 32
    100, 150, 5;                     \\# [10, 500] best 40
    1100, 1500, 100;                \\#[1000, 2000] best 540
    1100, 1300, 50;                 \\#[1000, 2000] best 240
    1800, 2000, 100;                 \\#[1000, 2000] best 620
    2000, 2800, 200;                 \\#[1000, 2000] best 750
    1800, 2400, 200;                 \\#[1000, 2000] best 250
    1000, 2000, 200
];

}
