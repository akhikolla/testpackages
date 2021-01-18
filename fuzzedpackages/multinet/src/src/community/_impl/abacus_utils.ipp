namespace uu {
namespace net {

template <typename M, typename L>
std::unique_ptr<CommunityStructure<PillarCommunity<L>>>
eclat_merge(
    const M* mnet,
    const std::unordered_map<const L*, CommunityStructure<Community<const Vertex*>>*>& single_layer_communities,
    int min_actors,
    int min_layers
)
{

    std::unordered_map<int, std::vector<std::string> > transactions;

    for (auto pair: single_layer_communities)
    {
        std::string layer_str = std::to_string(mnet->layers()->index_of(pair.first));
        int community = 0;

        for (auto c: *pair.second)
        {
            std::string community_str = std::to_string(community);

            for (auto node: *c)
            {
                transactions[mnet->actors()->index_of(node)].push_back(community_str + ":" + layer_str);
            }

            community++;
        }
    }

    FILE* f_inp = std::tmpfile();

    if (!f_inp)
    {
        throw core::FileNotFoundException("Cannot open input tmp file");
    }

    FILE* f_out = std::tmpfile();

    if (!f_out)
    {
        throw core::FileNotFoundException("Cannot open output tmp file");
    }

    FILE* f_tid = std::tmpfile();

    if (!f_tid)
    {
        throw core::FileNotFoundException("Cannot open output tmp file for transactions");
    }

    int     k = 0;             /* result variable */
    CCHAR   *recseps = NULL;      /* record  separators */
    CCHAR   *fldseps = NULL;      /* field   separators */
    CCHAR   *blanks  = NULL;      /* blank   characters */
    CCHAR   *comment = NULL;      /* comment characters */
    CCHAR   *hdr     = "";        /* record header  for output */
    CCHAR   *sep     = " ";       /* item separator for output */
    CCHAR   *imp     = " <- ";    /* implication sign for ass. rules */
    CCHAR   *dflt    = " (%S)";   /* default format for check */
    CCHAR   *info    = dflt;      /* format for information output */
    int     target   = ISR_CLOSED;       /* target type (e.g. closed/maximal) */
    ITEM    zmax     = ITEM_MAX;  /* maximum rule/item set size */
    double  smax     = 100;       /* maximum support of an item set */
    double  conf     = 80;        /* minimum confidence (in percent) */
    int     eval     = 'x';       /* additional evaluation measure */
    int     agg      = 'x';       /* aggregation mode for eval. measure */
    double  thresh   = 10;        /* threshold for evaluation measure */
    ITEM    prune    = ITEM_MIN;  /* (min. size for) evaluation pruning */
    int     sort     = 2;         /* flag for item sorting and recoding */
    int     algo     = 'a';       /* variant of eclat algorithm */
    int     mode     = ECL_DEFAULT|ECL_PREFMT;   /* search mode */
    int     pack     = 16;        /* number of bit-packed items */
    int     cmfilt   = -1;        /* mode for closed/maximal filtering */
    int     mtar     = 0;         /* mode for transaction reading */
    int     scan     = 0;         /* flag for scanable item output */
    ITEM    m;                    /* number of items */
    TID     n;                    /* number of transactions */

    TABREAD  *tread  = NULL; /* table/transaction reader */
    ITEMBASE *ibase  = NULL; /* item base */
    TABAG    *tabag  = NULL; /* transaction bag/multiset */
    ISREPORT *report = NULL; /* item set reporter */
    ECLAT    *eclat  = NULL; /* eclat miner object */

    for (size_t i=0; i<mnet->actors()->size(); i++)
    {
        for (std::string item: transactions[i])
        {
            fprintf(f_inp,"%s ",item.c_str());
        }

        fprintf(f_inp,"\n");
    }

    rewind(f_inp);

    if ((cmfilt >= 0) && (target & (ISR_CLOSED|ISR_MAXIMAL)))
    {
        mode |= (cmfilt > 0) ? ECL_VERT : ECL_HORZ;
    }

    mode |= ECL_TIDS;
    mode = (mode & ~ECL_FIM16)    /* add packed items to search mode */
           | ((pack <= 0) ? 0 : (pack < 16) ? pack : 16);
    mode |= ECL_VERBOSE|ECL_NOCLEAN;

    /* --- read item selection/appearance indicators --- */
    ibase = ib_create(0, 0);      /* create an item base */

    if (!ibase)
    {
        throw core::ExternalLibException("Cannot create item base");
    }

    tread = trd_create();         /* create a transaction reader */

    if (!tread)
    {
        throw core::ExternalLibException("Cannot create a transaction reader");
    }

    trd_allchs(tread, recseps, fldseps, blanks, "", comment);

    /* --- read transaction database --- */
    tabag = tbg_create(ibase);    /* create a transaction bag */

    if (!tabag)
    {
        throw core::ExternalLibException("Cannot create a transaction bag");
    }

    if (trd_open(tread, f_inp, NULL) != 0)
    {
        throw core::ExternalLibException("Cannot open the input transactions file ");
    }

    k = tbg_read(tabag, tread, mtar);

    if (k)
    {
        throw core::ExternalLibException("Cannot execute eternal function tbg_read");
    }

    trd_delete(tread, 1);         /* read the transaction database, */
    tread = NULL;                 /* then delete the table reader */
    m = ib_cnt(ibase);            /* get the number of items, */
    n = tbg_cnt(tabag);           /* the number of transactions, */

    if ((m <= 0) || (n <= 0))     /* check for at least one item */
    {
        return std::make_unique<CommunityStructure<PillarCommunity<L>>>();
    }

    /* --- find frequent item sets/association rules --- */
    eclat = eclat_create(target, -min_actors, smax, conf, min_layers, zmax,
                         eval, agg, thresh, algo, mode);

    if (!eclat)
    {
        throw core::ExternalLibException("Cannot create eclat miner");
    }

    k = eclat_data(eclat, tabag, 0, sort);

    if (k)
    {
        throw core::ExternalLibException("Cannot prepare data for eclat ");
    }

    report = isr_create(ibase);   /* create an item set reporter */

    if (!report)
    {
        throw core::ExternalLibException("Cannot configure report ");
    }

    k = eclat_report(eclat, report);

    if (k)
    {
        throw core::ExternalLibException("Cannot prepare report ");
    }

    if (isr_setfmt(report, scan, hdr, sep, imp, info) != 0)
    {
        throw core::ExternalLibException("Cannot set the oputput format string");
    }

    k = isr_tidopen(report, f_tid, NULL);  /* open the file for */

    if (k)
    {
        throw core::ExternalLibException("Cannot open transaction ids");
    }

    k = isr_open(report, f_out, NULL);

    if (k)
    {
        throw core::ExternalLibException("Cannot open report");
    }

    if (isr_setup(report))    /* open the item set file and */
    {
        throw core::ExternalLibException("Cannot setup itemset reporter");
    }

    k = eclat_mine(eclat, prune, 0);

    if (k)
    {
        throw core::ExternalLibException("Cannot run eclat_mine");
    }

    //if (stats)                    /* print item set statistics */
    //    isr_prstats(report, stdout, 0);
    isr_flush(report);
    isr_tidflush(report);

    rewind(report->file);
    rewind(report->tidfile);
    auto result = read_eclat_communities<M,L>(mnet, report->file, report->tidfile);

    /*if (f_inp!=NULL) {
     std::fclose(f_inp);
     f_inp = NULL;
     }
     if (f_out!=NULL) {
     std::fclose(f_out);
     f_out = NULL;
     }
     if (f_tid!=NULL) {
     std::fclose(f_tid);
     f_tid=NULL;
     }*/

    if (isr_close   (report) != 0)
    {
        throw core::ExternalLibException("Cannot close itemset output file");
    }

    if (isr_tidclose(report) != 0)
    {
        throw core::ExternalLibException("Cannot close transaction id output file");
    }

    //CLEANUP;
    if (eclat)
    {
        eclat_delete(eclat, 0);
    }

    if (report)
    {
        isr_delete(report, 0);
    }

    if (tabag)
    {
        tbg_delete(tabag,  0);
    }

    if (tread)
    {
        trd_delete(tread,  1);
    }

    if (ibase)
    {
        ib_delete (ibase);
    }

    return result;
}


template <typename M, typename L>
void
read_layers(
    const M* mnet,
    PillarCommunity<L>* com,
    FILE* file
)
{
    bool found_separator=false;
    bool reading_number=false;
    int c;
    int lid; // layer identifier

    while (1)
    {
        c=getc(file);

        if (c=='\n' || c==EOF)
        {
            if (reading_number)
            {
                com->add_layer(mnet->layers()->at(lid));
            }

            return;
        }

        if (!found_separator && c!=':')
        {
            continue;
        }

        else if (c==':')
        {
            found_separator=true;
            continue;
        }

        if (c>='0' && c<='9')
        {
            if (!reading_number)
            {
                lid=(int)(c-'0');
            }

            else
            {
                lid=lid*10+(int)(c-'0');
            }

            reading_number=true;
        }

        else if (c==' ')
        {
            if (reading_number)
            {
                com->add_layer(mnet->layers()->at(lid));
            }

            reading_number=false;
            found_separator=false;
        }
    }
}

template <typename M, typename L>
int
read_actors(
    const M* mnet,
    PillarCommunity<L>* com,
    FILE* tidfile
)
{
    size_t id;
    int c;
    bool reading_number=false;

    while (1)
    {
        c=getc(tidfile);

        if (c>='0' && c<='9')
        {
            if (!reading_number)
            {
                id=(int)(c-'0');
            }

            else
            {
                id=id*10+(int)(c-'0');
            }

            reading_number=true;
        }

        else if (c==' ' || c=='\n' || c==EOF)
        {
            if (reading_number)
            {
                com->add_actor(mnet->actors()->at(id-1));
            }

            reading_number=false;

            if (c==EOF)
            {
                return 0;
            }

            if (c=='\n')
            {
                return 1;
            }
        }
    }

}

template <typename M, typename L>
std::unique_ptr<CommunityStructure<PillarCommunity<L>>>
read_eclat_communities(
    const M* mnet,
    FILE* file,
    FILE* tidfile
)
{
    auto communities = std::make_unique<CommunityStructure<PillarCommunity<L>>>();
    auto current = std::make_unique<PillarCommunity<L>>();

    while (read_actors(mnet,current.get(),tidfile))
    {
        read_layers(mnet,current.get(),file);

        if (current->get_actors().size()>0)
        {
            communities->add(std::move(current));
            current = std::make_unique<PillarCommunity<L>>();
        }
    }

    communities->add(std::move(current));
    return communities;
}

}
}

