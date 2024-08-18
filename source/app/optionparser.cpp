/*
 * The Lean Mean C++ Option Parser
 *
 * Copyright (C) 2012 Matthias S. Benkmann
 *
 * The "Software" in the following 2 paragraphs refers to this file containing
 * the code to The Lean Mean C++ Option Parser.
 * The "Software" does NOT refer to any other files which you
 * may have received alongside this file (e.g. as part of a larger project that
 * incorporates The Lean Mean C++ Option Parser).
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software, to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to permit
 * persons to whom the Software is furnished to do so, subject to the following
 * conditions:
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "app/optionparser/optionparser.h"

namespace option
{

int Option::type() const
{
    return desc == 0 ? 0 : desc->type;
}

int Option::index() const
{
    return desc == 0 ? -1 : (int)desc->index;
}

int Option::count()
{
    int c = ( desc == 0 ? 0 : 1 );
    Option *p = first();
    while( !p->isLast() )
    {
        ++c;
        p = p->next_;
    };
    return c;
}

bool Option::isFirst() const
{
    return isTagged( prev_ );
}

bool Option::isLast() const
{
    return isTagged( next_ );
}

Option *Option::first()
{
    Option *p = this;
    while( !p->isFirst() )
        p = p->prev_;
    return p;
}

Option *Option::last()
{
    return first()->prevwrap();
}

Option *Option::prev()
{
    return isFirst() ? 0 : prev_;
}

Option *Option::prevwrap()
{
    return untag( prev_ );
}

Option *Option::next()
{
    return isLast() ? 0 : next_;
}

const Option *Option::next() const
{
    return isLast() ? 0 : next_;
}

Option *Option::nextwrap()
{
    return untag( next_ );
}

void Option::append( Option *new_last )
{
    Option *p = last();
    Option *f = first();
    p->next_ = new_last;
    new_last->prev_ = p;
    new_last->next_ = tag( f );
    f->prev_ = tag( new_last );
}

Option::Option() : desc( 0 ), name( 0 ), arg( 0 ), namelen( 0 )
{
    prev_ = tag( this );
    next_ = tag( this );
}

Option::Option( const Descriptor *desc_, const char *name_, const char *arg_ )
{
    init( desc_, name_, arg_ );
}

void Option::operator=( const Option &orig )
{
    init( orig.desc, orig.name, orig.arg );
}

Option::Option( const Option &orig )
{
    init( orig.desc, orig.name, orig.arg );
}

void Option::init( const Descriptor *desc_, const char *name_, const char *arg_ )
{
    desc = desc_;
    name = name_;
    arg = arg_;
    prev_ = tag( this );
    next_ = tag( this );
    namelen = 0;
    if( name == 0 )
        return;
    namelen = 1;
    if( name[0] != '-' )
        return;
    while( name[namelen] != 0 && name[namelen] != '=' )
        ++namelen;
}

Option *Option::tag( Option *ptr )
{
    return (Option *)( (unsigned long long)ptr | 1 );
}

Option *Option::untag( Option *ptr )
{
    return (Option *)( (unsigned long long)ptr & ~1ull );
}

bool Option::isTagged( Option *ptr )
{
    return ( (unsigned long long)ptr & 1 );
}

Stats::Stats() : buffer_max( 1 ), options_max( 1 ) // 1 more than necessary as sentinel
{
}

Stats::Stats(
    bool gnu, const Descriptor usage[], int argc, const char **argv, int min_abbr_len, bool single_minus_longopt )
    : buffer_max( 1 ), options_max( 1 ) // 1 more than necessary as sentinel
{
    add( gnu, usage, argc, argv, min_abbr_len, single_minus_longopt );
}

Stats::Stats( bool gnu, const Descriptor usage[], int argc, char **argv, int min_abbr_len, bool single_minus_longopt )
    : buffer_max( 1 ), options_max( 1 ) // 1 more than necessary as sentinel
{
    add( gnu, usage, argc, (const char **)argv, min_abbr_len, single_minus_longopt );
}

Stats::Stats( const Descriptor usage[], int argc, const char **argv, int min_abbr_len, bool single_minus_longopt )
    : buffer_max( 1 ), options_max( 1 ) // 1 more than necessary as sentinel
{
    add( false, usage, argc, argv, min_abbr_len, single_minus_longopt );
}

Stats::Stats( const Descriptor usage[], int argc, char **argv, int min_abbr_len, bool single_minus_longopt )
    : buffer_max( 1 ), options_max( 1 ) // 1 more than necessary as sentinel
{
    add( false, usage, argc, (const char **)argv, min_abbr_len, single_minus_longopt );
}

void Stats::add(
    bool gnu, const Descriptor usage[], int argc, const char **argv, int min_abbr_len, bool single_minus_longopt )
{
    // determine size of options array. This is the greatest index used in the usage + 1
    int i = 0;
    while( usage[i].shortopt != 0 )
    {
        if( usage[i].index + 1 >= options_max )
            options_max = ( usage[i].index + 1 ) + 1; // 1 more than necessary as sentinel

        ++i;
    }

    CountOptionsAction action( &buffer_max );
    Parser::workhorse( gnu, usage, argc, argv, action, single_minus_longopt, false, min_abbr_len );
}

void Stats::add(
    bool gnu, const Descriptor usage[], int argc, char **argv, int min_abbr_len, bool single_minus_longopt )
{
    add( gnu, usage, argc, (const char **)argv, min_abbr_len, single_minus_longopt );
}

void Stats::add( const Descriptor usage[], int argc, const char **argv, int min_abbr_len, bool single_minus_longopt )
{
    add( false, usage, argc, argv, min_abbr_len, single_minus_longopt );
}

void Stats::add( const Descriptor usage[], int argc, char **argv, int min_abbr_len, bool single_minus_longopt )
{
    add( false, usage, argc, (const char **)argv, min_abbr_len, single_minus_longopt );
}

void Parser::parse( bool gnu,
                    const Descriptor usage[],
                    int argc,
                    const char **argv,
                    Option options[],
                    Option buffer[],
                    int min_abbr_len,
                    bool single_minus_longopt,
                    int bufmax )
{
    StoreOptionAction action( *this, options, buffer, bufmax );
    err = !workhorse( gnu, usage, argc, argv, action, single_minus_longopt, true, min_abbr_len );
}

bool Parser::workhorse( bool gnu,
                        const Descriptor usage[],
                        int numargs,
                        const char **args,
                        Action &action,
                        bool single_minus_longopt,
                        bool print_errors,
                        int min_abbr_len )
{
    // protect against NULL pointer
    if( args == 0 )
        numargs = 0;

    int nonops = 0;

    while( numargs != 0 && *args != 0 )
    {
        const char *param = *args; // param can be --long-option, -srto or non-option argument

        // in POSIX mode the first non-option argument terminates the option list
        // a lone minus character is a non-option argument
        if( param[0] != '-' || param[1] == 0 )
        {
            if( gnu )
            {
                ++nonops;
                ++args;
                if( numargs > 0 )
                    --numargs;
                continue;
            }
            else
                break;
        }

        // -- terminates the option list. The -- itself is skipped.
        if( param[1] == '-' && param[2] == 0 )
        {
            shift( args, nonops );
            ++args;
            if( numargs > 0 )
                --numargs;
            break;
        }

        bool handle_short_options;
        const char *longopt_name;
        if( param[1] == '-' ) // if --long-option
        {
            handle_short_options = false;
            longopt_name = param + 2;
        }
        else
        {
            handle_short_options = true;
            longopt_name = param + 1; // for testing a potential -long-option
        }

        bool try_single_minus_longopt = single_minus_longopt;
        bool have_more_args = ( numargs > 1 || numargs < 0 ); // is referencing argv[1] valid?

        do // loop over short options in group, for long options the body is executed only once
        {
            int idx = 0;

            const char *optarg = 0;

            /******************** long option **********************/
            if( handle_short_options == false || try_single_minus_longopt )
            {
                idx = 0;
                while( usage[idx].longopt != 0 && !streq( usage[idx].longopt, longopt_name ) )
                    ++idx;

                if( usage[idx].longopt == 0 && min_abbr_len > 0 ) // if we should try to match abbreviated long options
                {
                    int i1 = 0;
                    while( usage[i1].longopt != 0 && !streqabbr( usage[i1].longopt, longopt_name, min_abbr_len ) )
                        ++i1;
                    if( usage[i1].longopt != 0 )
                    { // now test if the match is unambiguous by checking for another match
                        int i2 = i1 + 1;
                        while( usage[i2].longopt != 0 && !streqabbr( usage[i2].longopt, longopt_name, min_abbr_len ) )
                            ++i2;

                        if( usage[i2].longopt ==
                            0 ) // if there was no second match it's unambiguous, so accept i1 as idx
                            idx = i1;
                    }
                }

                // if we found something, disable handle_short_options (only relevant if single_minus_longopt)
                if( usage[idx].longopt != 0 )
                    handle_short_options = false;

                try_single_minus_longopt = false; // prevent looking for longopt in the middle of shortopt group

                optarg = longopt_name;
                while( *optarg != 0 && *optarg != '=' )
                    ++optarg;
                if( *optarg == '=' ) // attached argument
                    ++optarg;
                else
                    // possibly detached argument
                    optarg = ( have_more_args ? args[1] : 0 );
            }

            /************************ short option ***********************************/
            if( handle_short_options )
            {
                if( *++param == 0 ) // point at the 1st/next option character
                    break;          // end of short option group

                idx = 0;
                while( usage[idx].shortopt != 0 && !instr( *param, usage[idx].shortopt ) )
                    ++idx;

                if( param[1] == 0 ) // if the potential argument is separate
                    optarg = ( have_more_args ? args[1] : 0 );
                else
                    // if the potential argument is attached
                    optarg = param + 1;
            }

            const Descriptor *descriptor = &usage[idx];

            if( descriptor->shortopt == 0 ) /**************  unknown option ********************/
            {
                // look for dummy entry (shortopt == "" and longopt == "") to use as Descriptor for unknown options
                idx = 0;
                while( usage[idx].shortopt != 0 && ( usage[idx].shortopt[0] != 0 || usage[idx].longopt[0] != 0 ) )
                    ++idx;
                descriptor = ( usage[idx].shortopt == 0 ? 0 : &usage[idx] );
            }

            if( descriptor != 0 )
            {
                Option option( descriptor, param, optarg );
                switch( descriptor->check_arg( option, print_errors ) )
                {
                case ARG_ILLEGAL:
                    return false; // fatal
                case ARG_OK:
                    // skip one element of the argument vector, if it's a separated argument
                    if( optarg != 0 && have_more_args && optarg == args[1] )
                    {
                        shift( args, nonops );
                        if( numargs > 0 )
                            --numargs;
                        ++args;
                    }

                    // No further short options are possible after an argument
                    handle_short_options = false;

                    break;
                case ARG_IGNORE:
                case ARG_NONE:
                    option.arg = 0;
                    break;
                }

                if( !action.perform( option ) )
                    return false;
            }

        } while( handle_short_options );

        shift( args, nonops );
        ++args;
        if( numargs > 0 )
            --numargs;

    } // while

    if( numargs > 0 && *args == 0 ) // It's a bug in the caller if numargs is greater than the actual number
        numargs = 0;                // of arguments, but as a service to the user we fix this if we spot it.

    if( numargs < 0 ) // if we don't know the number of remaining non-option arguments
    {                 // we need to count them
        numargs = 0;
        while( args[numargs] != 0 )
            ++numargs;
    }

    return action.finished( numargs + nonops, args - nonops );
}
} // namespace option
// namespace option
