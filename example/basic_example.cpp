// internal
#include <promille/promille.hpp>

#include "../test/source/dummy_residual_model.hpp"

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
    int verbose_flag {0};
    int print_flag {0};

    option long_options[] = {/* These options set a flag. */
                             {"verbose", no_argument, &verbose_flag, 1},
                             {"brief", no_argument, &verbose_flag, 0},
                             {"print", no_argument, 0, 'p'},
                             {nullptr, 0, nullptr, 0}};

    while (true) {
        /* getopt_long stores the option index here. */
        int option_index = 0;

        const int c = getopt_long(argc, argv, "p", long_options, &option_index);

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

            case 'p':
                print_flag = 1;
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

    promille::promille mille("test_", "test.bin");

    mille.set_verbose(2);

    // FIRST WAY: define parameters and pass their indices to add_plane()
    // Useful, if these params must be shared between different measurement planes.
    mille.add_global_parameter(11, -40, "Ra1");
    mille.add_global_parameter(12, -50, "Rb1");
    mille.add_global_parameter(13, -60, "Rc1");

    auto straw_planes1 = mille.make_model_planes<promille_tests::dummy_residual_model<float, float>>();

    straw_planes1.add_plane(1, 11, 12, 13)
        .set_globals_configuration(promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FIXED)
        .set_locals_configuration(promille::Kind::FREE, promille::Kind::FIXED);

    // SECOND WAY: define parameters inline, useful if not needed anymore
    straw_planes1
        .add_plane(2,
                   mille.add_global_parameter(21, -40, "Ra2"),
                   mille.add_global_parameter(22, -50, "Rb2"),
                   mille.add_global_parameter(23, -60, "Rc2"))
        .set_globals_configuration(promille::Kind::FIXED, promille::Kind::FREE, promille::Kind::FIXED)
        .set_locals_configuration(promille::Kind::FREE, promille::Kind::FIXED);

    std::array<promille::Kind, 3> plane_3_global_kinds = {promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FREE};
    std::array<promille::Kind, 2> plane_3_local_kinds = {
        promille::Kind::FREE,
        promille::Kind::FREE,
    };

    straw_planes1
        .add_plane(3,
                   mille.add_global_parameter(31, -40, "Ra3"),
                   mille.add_global_parameter(32, -50, "Rb3"),
                   mille.add_global_parameter(33, -60, "Rc3"))
        .set_globals_configuration(plane_3_global_kinds)
        .set_locals_configuration(plane_3_local_kinds);

    // FOURTH WAY: mix different ways to initialize
    mille.add_global_parameter(41, -40, "Ra4");
    mille.add_global_parameter(42, -50, "Rb4");
    mille.add_global_parameter(43, -60, "Rc4");

    straw_planes1.add_plane(4, 41, 42, 43)
        .set_globals_configuration(promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FREE)
        .set_locals_configuration(promille::Kind::FREE, promille::Kind::FIXED);

    if (print_flag) {
        mille.print();
        straw_planes1.print();
    }

    auto meas = straw_planes1.plane(1).add_measurement(0);

    meas.model().print();

    return 0;
}
