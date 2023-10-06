// internal
#include <promille/promille.hpp>

#include "../test/source/dummy_alignment_model.hpp"

// ROOT
#include <Math/Point3D.h>
#include <Math/Vector3D.h>

// system
#include <string>

#include <getopt.h>

using ROOT::Math::XYZPoint;
using ROOT::Math::XYZVector;

auto main(int argc, char* argv[]) -> int
{
    /* Flag set by ‘--verbose’. */
    int verbose_flag;

    int c;

    while (1) {
        static struct option long_options[] = {/* These options set a flag. */
                                               {"verbose", no_argument, &verbose_flag, 1},
                                               {"brief", no_argument, &verbose_flag, 0},
                                               /* These options don’t set a flag.
                                                  We distinguish them by their indices. */
                                               // {"add", no_argument, 0, 'a'},
                                               // {"append", no_argument, 0, 'b'},
                                               // {"delete", required_argument, 0, 'd'},
                                               // {"create", required_argument, 0, 'c'},
                                               // {"file", required_argument, 0, 'f'},
                                               {0, 0, 0, 0}};
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long(argc, argv, "abc:d:f:", long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c) {
            case 0:
                /* If this option set a flag, do nothing else now. */
                if (long_options[option_index].flag != 0)
                    break;
                printf("option %s", long_options[option_index].name);
                if (optarg)
                    printf(" with arg %s", optarg);
                printf("\n");
                break;

            case 'a':
                puts("option -a\n");
                break;

            case 'b':
                puts("option -b\n");
                break;

            case 'c':
                printf("option -c with value `%s'\n", optarg);
                break;

            case 'd':
                printf("option -d with value `%s'\n", optarg);
                break;

            case 'f':
                printf("option -f with value `%s'\n", optarg);
                break;

            case '?':
                /* getopt_long already printed an error message. */
                break;

            default:
                abort();
        }
    }

    /* Instead of reporting ‘--verbose’
       and ‘--brief’ as they are encountered,
       we report the final status resulting from them. */
    if (verbose_flag)
        puts("verbose flag is set");

    std::string input_data_file;

    /* Print any remaining command line arguments (not options). */
    if (optind < argc) {
        printf("non-option ARGV-elements: ");
        while (optind < argc) {
            printf("%s ", argv[optind]);
            input_data_file = std::string(argv[optind++]);
            break;
        }
        putchar('\n');
    }

    promille::promille<mb_tests::dummy_residual_model<float, promille::euler::zyz>> mille("test_", "test");

    mille.set_local_parameter(0, "Bx");
    mille.set_local_parameter(1, "By");
    mille.set_local_parameter(2, "Tx");
    mille.set_local_parameter(3, "Ty");

    mille.add_global_parameter(11, -40, "Ra1");
    mille.add_global_parameter(12, -50, "Rb1");
    mille.add_global_parameter(13, -60, "Rc1");
    mille.add_global_parameter(14, -30, "Tx1");
    mille.add_global_parameter(15, -20, "Ty1");
    mille.add_global_parameter(16, -10, "Tz1");

    mille.add_global_parameter(21, -40, "Ra2");
    mille.add_global_parameter(22, -50, "Rb2");
    mille.add_global_parameter(23, -60, "Rc2");
    mille.add_global_parameter(24, -30, "Tx2");
    mille.add_global_parameter(25, -20, "Ty2");
    mille.add_global_parameter(26, -10, "Tz2");

    mille.add_global_parameter(31, -40, "Ra3");
    mille.add_global_parameter(32, -50, "Rb3");
    mille.add_global_parameter(33, -60, "Rc3");
    mille.add_global_parameter(34, -30, "Tx3");
    mille.add_global_parameter(35, -20, "Ty3");
    mille.add_global_parameter(36, -10, "Tz3");

    mille.add_global_parameter(41, -40, "Ra4");
    mille.add_global_parameter(42, -50, "Rb4");
    mille.add_global_parameter(43, -60, "Rc4");
    mille.add_global_parameter(44, -30, "Tx4");
    mille.add_global_parameter(45, -20, "Ty4");
    mille.add_global_parameter(46, -10, "Tz4");

    mille.add_plane(1, 11, 12, 13, 14, 15, 16)
        .set_globals_configuration(promille::Kind::FREE,
                                   promille::Kind::FREE,
                                   promille::Kind::FIXED,
                                   promille::Kind::FIXED,
                                   promille::Kind::FIXED,
                                   promille::Kind::FREE)
        .set_locals_configuration(promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FREE);

    mille.add_plane(2, 21, 22, 23, 24, 25, 26)
        .set_globals_configuration(promille::Kind::FREE,
                                   promille::Kind::FREE,
                                   promille::Kind::FIXED,
                                   promille::Kind::FIXED,
                                   promille::Kind::FIXED,
                                   promille::Kind::FREE)
        .set_locals_configuration(promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FREE);

    mille.add_plane(3, 31, 32, 33, 34, 35, 36)
        .set_globals_configuration(promille::Kind::FREE,
                                   promille::Kind::FREE,
                                   promille::Kind::FIXED,
                                   promille::Kind::FIXED,
                                   promille::Kind::FIXED,
                                   promille::Kind::FREE)
        .set_locals_configuration(promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FREE);

    mille.add_plane(4, 41, 42, 43, 44, 45, 46)
        .set_globals_configuration(promille::Kind::FREE,
                                   promille::Kind::FREE,
                                   promille::Kind::FIXED,
                                   promille::Kind::FIXED,
                                   promille::Kind::FIXED,
                                   promille::Kind::FREE)
        .set_locals_configuration(promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FREE);

    mille.print(true);

    mille.add_measurement(1, 0, XYZPoint(0, 0, 1), XYZVector(2, 3, 4));

    return 0;
}
