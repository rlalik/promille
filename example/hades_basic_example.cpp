#include <string>

#include <StrawAlignment/StrawAlignment.hpp>
#include <TMath.h>
#include <getopt.h>

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

    auto mille = SA::MilleBuilder<SA::euler::zyz>("test_", "hades_basic_example_alignment_data.bin");
    mille.set_verbose(verbose_flag);

    if (input_data_file.length()) {
        mille.read_global_data(input_data_file);
    } else {
        mille.add_planes_globals({0, SA::Kind::FREE},
                                 {0, SA::Kind::FIXED},
                                 {0, SA::Kind::FIXED},
                                 {TMath::DegToRad() * 0, SA::Kind::FIXED},
                                 {0, SA::Kind::FIXED},
                                 {0, SA::Kind::FIXED}, 0, 0, 0, 0, 0, 0);

        mille.add_planes_globals({0, SA::Kind::FIXED},
                                 {0, SA::Kind::FREE},
                                 {0, SA::Kind::FIXED},
                                 {TMath::DegToRad() * 90, SA::Kind::FIXED},
                                 {0, SA::Kind::FIXED},
                                 {0, SA::Kind::FIXED}, 0, 0, 0, 0, 0, 0);

        mille.add_planes_globals({0, SA::Kind::FIXED},
                                 {0, SA::Kind::FREE},
                                 {0, SA::Kind::FIXED},
                                 {TMath::DegToRad() * 90, SA::Kind::FIXED},
                                 {0, SA::Kind::FIXED},
                                 {0, SA::Kind::FIXED}, 0, 0, 0, 0, 0, 0);

        mille.add_planes_globals({0, SA::Kind::FREE},
                                 {0, SA::Kind::FIXED},
                                 {0, SA::Kind::FIXED},
                                 {TMath::DegToRad() * 0, SA::Kind::FIXED},
                                 {0, SA::Kind::FIXED},
                                 {0, SA::Kind::FIXED}, 0, 0, 0, 0, 0, 0);

        mille.add_planes_globals({0, SA::Kind::FREE},
                                 {0, SA::Kind::FIXED},
                                 {0, SA::Kind::FIXED},
                                 {TMath::DegToRad() * 0, SA::Kind::FIXED},
                                 {0, SA::Kind::FIXED},
                                 {0, SA::Kind::FIXED}, 0, 0, 0, 0, 0, 0);

        mille.add_planes_globals({0, SA::Kind::FIXED},
                                 {0, SA::Kind::FREE},
                                 {0, SA::Kind::FIXED},
                                 {TMath::DegToRad() * 90, SA::Kind::FIXED},
                                 {0, SA::Kind::FIXED},
                                 {0, SA::Kind::FIXED}, 0, 0, 0, 0, 0, 0);

        mille.add_planes_globals({0, SA::Kind::FREE},
                                 {0, SA::Kind::FREE},
                                 {0, SA::Kind::FIXED},
                                 {TMath::DegToRad() * 45, SA::Kind::FIXED},
                                 {0, SA::Kind::FIXED},
                                 {0, SA::Kind::FIXED}, 0, 0, 0, 0, 0, 0);

        mille.add_planes_globals({0, SA::Kind::FREE},
                                 {0, SA::Kind::FREE},
                                 {0, SA::Kind::FIXED},
                                 {TMath::DegToRad() * -45, SA::Kind::FIXED},
                                 {0, SA::Kind::FIXED},
                                 {0, SA::Kind::FIXED}, 0, 0, 0, 0, 0, 0);
    }

    mille.print(true);

    mille.add_local(0, 0, 1, 2, 3, 4, 0, 0, 0, 0.1, 0.15);
    mille.add_local(1, 1, 2, 3, 4, 5, 0, 0, 0, 0.1, 0.15);
    mille.add_local(2, 2, 3, 4, 5, 6, 0, 0, 0, 0.1, 0.15);
    mille.add_local(3, 0, 1, 2, 3, 4, 0, 0, 0, 0.1, 0.15);
    mille.add_local(4, 1, 2, 3, 4, 5, 0, 0, 0, 0.1, 0.15);
    mille.add_local(5, 2, 3, 4, 5, 6, 0, 0, 0, 0.1, 0.15);
    mille.end();

    mille.add_local(0, 1, 2, 3, 4, 5, 0, 0, 0, 0.1, 0.15);
    mille.add_local(1, 2, 3, 4, 5, 6, 0, 0, 0, 0.1, 0.15);
    mille.add_local(2, 3, 4, 5, 6, 7, 0, 0, 0, 0.1, 0.15);
    mille.end();

    return 0;
}
